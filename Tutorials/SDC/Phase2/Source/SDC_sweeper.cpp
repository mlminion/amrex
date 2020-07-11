#include "myfunc.H"
#include "myfunc_F.H"
#include <AMReX_BCUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <iostream>
#include <AMReX_PlotFileUtil.H>

void SDC_advance(MultiFab& phi_old,
		 MultiFab& phi_new,
		 std::array<MultiFab, AMREX_SPACEDIM>& flux,
		 Real dt,
		 const Geometry& geom,
		 const Vector<BCRec>& bc,
		 MLMG&  mlmg,
		 MLABecLaplacian& mlabec,
		 SDCstruct &SDC, 
		 std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
		 std::array<MultiFab,AMREX_SPACEDIM>& prod_stor,
		 Real time, 
		 MultiFab& bdry_values,
		 int &tot_Vcycle,int &tot_res_iter,int &tot_SDC_sweep)
{

  /*  This is a multi-implicit SDC example time step for an 
      advection-diffusion-reaction equation of the form
      
      phi_t = A(phi)+D(phi)+R(phi)
      
      The advection is treated explicilty and the diffusion and reaction implicitly
      and uncoupled
      
      the constants a,d, and r control the strength of each term  */ 
    
  Real qij;   //  Temp used for quadrature weight
  Real current_time = time;  //  Updated at SDC nodes in loop
  ParmParse pp;
  int Nprob;
  Real k_freq;
  pp.get("Nprob",Nprob);
  pp.get("k_freq",k_freq);
  // Copy old phi into first SDC node
  MultiFab::Copy(SDC.sol[0],phi_old, 0, 0, 1, phi_old.nGrow());
  
  // Fill the ghost cells of each grid from the other grids
  // includes periodic domain boundaries
  // Fill Dirichlet values
  bdry_values.setVal(0,0); bdry_values.setBndry(0);
  if(Nprob<3){
    for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
      {
	const Box& bx = mfi.validbox();
	fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
			 BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
			 geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&current_time, &Nprob);
      }
    
    mlabec.fourthOrderBCFill(SDC.sol[0],bdry_values);
  }
  // Fill periodic values and interal BC ghost cells
  SDC.sol[0].FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries (doesn't do Dirichlet; probably only one ghost cell).
  // FillDomainBoundary(SDC.sol[0], geom, bc);
  
  //  Compute the first function value (Need boundary filled here at current time explicit).
  int sdc_m=0;
  SDC_feval(flux,geom,bc,SDC,face_bcoef,prod_stor,sdc_m,-1,time);
  
  // Copy first function value to all nodes
  for (int sdc_n = 1; sdc_n < SDC.Nnodes; sdc_n++)
    {
      current_time = time+dt*SDC.qnodes[sdc_n];
      MultiFab::Copy(SDC.sol[sdc_n],SDC.sol[0], 0, 0, 1, SDC.sol[0].nGrow());
      // Subsequent Calculation doesn't need BC conditions ... fills out forcing term.  For time dependent BC, this would be inconsistent
      SDC_feval(flux,geom,bc,SDC,face_bcoef,prod_stor,sdc_n,0,current_time);
      
      //  Copy diffusion and reaction terms
      MultiFab::Copy(SDC.f[1][sdc_n],SDC.f[1][0], 0, 0, 1,SDC.f[1][0].nGrow() );
      if (SDC.Npieces==3)
	MultiFab::Copy(SDC.f[2][sdc_n],SDC.f[2][0], 0, 0, 1, SDC.f[2][0].nGrow());
    }


  //  Now do the actual sweeps
  int k=0;
  Real sdc_res=10.0e10;
  Real tol_abs_SDC=0.0;
  pp.query("tol_abs_SDC",tol_abs_SDC);  
  while ( (sdc_res > tol_abs_SDC) & (k <= SDC.Nsweeps) )  //  Loop over residual solves
    {
      amrex::Print() << "SDC sweep " << k << ", substep " << sdc_m+1 <<"---\n";
      //  Compute RHS integrals
      SDC.SDC_rhs_integrals(dt);
      
      //  Substep over SDC nodes
      for (sdc_m = 0; sdc_m < SDC.Nnodes-1; sdc_m++)
	{
	  amrex::Print() << "--- SDC substep " << sdc_m+1 <<"\n";
	  
	  // use phi_new as rhs and fill (overwrite) it with terms at this iteration
	  SDC.SDC_rhs_k_plus_one(phi_new,dt,sdc_m);
	  
	  // get the best initial guess for implicit solve
	  MultiFab::Copy(SDC.sol[sdc_m+1],phi_new, 0, 0, 1, phi_new.nGrow());
	  for ( MFIter mfi(SDC.sol[sdc_m+1]); mfi.isValid(); ++mfi )
	    {
	      //	      const Box& bx = mfi.validbox();
	      qij = dt*SDC.Qimp[sdc_m][sdc_m+1];
	      SDC.sol[sdc_m+1][mfi].saxpy(qij,SDC.f[1][sdc_m+1][mfi]);
	    }
	  // Get time value at next SDC node for implicit solve
	  current_time = time+dt*SDC.qnodes[sdc_m+1];
	  
	  // Fill Dirichlet Values
	  if (Nprob<3){
            for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
	      {
		const Box& bx = mfi.validbox();  //  Is 'const necessary here?
                fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
                                 geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&current_time, &Nprob);
	      }
            mlabec.fourthOrderBCFill(SDC.sol[sdc_m+1],bdry_values);
	  }
	  // Fill periodic values and interal ghost cells
	  SDC.sol[sdc_m+1].FillBoundary(geom.periodicity());
        
	  // Solve for the first implicit part
	  SDC_fcomp(phi_new, flux, geom, bc, SDC, mlmg, mlabec,dt,face_bcoef,prod_stor,sdc_m+1,1, current_time, tot_Vcycle, tot_res_iter);
	  
	  if (SDC.Npieces==3)
	    {
	      // Build rhs for 2nd solve
	      MultiFab::Copy(phi_new, SDC.sol[sdc_m+1],0, 0, 1, SDC.sol[sdc_m+1].nGrow());
	      
	      // Add in the part for the 2nd implicit term to rhs
	      SDC.SDC_rhs_misdc(phi_new,dt,sdc_m);
	      
	      // Solve for the second implicit part
	      SDC_fcomp(phi_new, flux, geom, bc, SDC, mlmg, mlabec,dt,face_bcoef,prod_stor,sdc_m+1,2, current_time,tot_Vcycle,tot_res_iter );
	    }
	  // Compute the function values at node sdc_m+1
	  
	  // Shouldn't need the following BC code as the solve should linearly maintain it.
	  if (Nprob<3){
	    for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
	      {          const Box& bx = mfi.validbox();
		fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
				 BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
				 geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&current_time, &Nprob);
	      }
	    mlabec.fourthOrderBCFill(SDC.sol[sdc_m+1],bdry_values);
	  }
	  SDC.sol[sdc_m+1].FillBoundary(geom.periodicity());
	  //       mlabec.fillSolutionBC(0, SDC.sol[sdc_m+1], &bdry_values);

	  //  Evaluate all the flux terms with new solution  
	  SDC_feval(flux,geom,bc,SDC,face_bcoef,prod_stor,sdc_m+1,-1,current_time);
	  
	} // end SDC substep loop (sdc_m)
      
      k++;
    }  // end sweeps loop (k)
  
  tot_SDC_sweep+=k-1;  //  increment total sweeps done
  // Return the last node in SDC.sol
  MultiFab::Copy(phi_new, SDC.sol[SDC.Nnodes-1], 0, 0, 1, SDC.sol[SDC.Nnodes-1].nGrow());

}

void SDC_feval(std::array<MultiFab, AMREX_SPACEDIM>& flux,
	       const Geometry& geom,
	       const Vector<BCRec>& bc,
	       SDCstruct &SDC,
	       std::array<MultiFab, AMREX_SPACEDIM>& face_bcoef,
	       std::array<MultiFab, AMREX_SPACEDIM>& prod_stor,
	       int sdc_m,int npiece, Real time)
{
  /*  Evaluate explicitly the rhs terms of the equation at the SDC node "sdc_m".
      The input parameter "npiece" describes which term to do.  
      If npiece = -1, do all the pieces */
  ParmParse pp;
  const Box& domain_bx = geom.Domain();
  const Real* dx = geom.CellSize();
  int nlo,nhi;

  //  Read in the coefficients for A-D-R
  Real a,d,r;
  pp.query("a",a);
  pp.query("d",d);
  pp.query("r",r);
  int Nprob,Lord;
  Real k_freq;
  Real epsilon=0.25;
  pp.get("Nprob",Nprob);
  pp.get("Lord",Lord);
  pp.get("k_freq",k_freq);
  pp.get("epsilon",epsilon);

  // Decide which pieces to loop over
  if (npiece < 0)
    {
      nlo=0;
      nhi=SDC.Npieces;
    }
  else
    {
      nlo=npiece;
      nhi=npiece+1;
    }
    
  for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();
      for (int n = nlo; n < nhi; n++)    
	{
	  SDC_feval_F(BL_TO_FORTRAN_BOX(bx),
		      BL_TO_FORTRAN_BOX(domain_bx),
		      BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
		      BL_TO_FORTRAN_ANYD(flux[0][mfi]),
		      BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
		      BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
		      BL_TO_FORTRAN_ANYD(SDC.f[n][sdc_m][mfi]),
		      dx,&a,&d,&r,
		      BL_TO_FORTRAN_ANYD(face_bcoef[0][mfi]),
		      BL_TO_FORTRAN_ANYD(face_bcoef[1][mfi]),
		      BL_TO_FORTRAN_ANYD(prod_stor[0][mfi]),
		      BL_TO_FORTRAN_ANYD(prod_stor[1][mfi]),
		      &n, &time, &epsilon, &k_freq, &Nprob, &Lord);
	}
      
    }
}
void SDC_fcomp(MultiFab& rhs,
	       std::array<MultiFab, AMREX_SPACEDIM>& flux,
	       const Geometry& geom,
	       const Vector<BCRec>& bc,
	       SDCstruct &SDC,
	       MLMG &mlmg,
	       MLABecLaplacian& mlabec,	      
	       Real dt,
	       std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
	       std::array<MultiFab,AMREX_SPACEDIM>& prod_stor,
	       int sdc_m,int npiece, Real time,
	       int &tot_Vcycle, int &tot_res_iter)
{
  /*  Solve implicitly for the  terms of the equation at the SDC node "sdc_m".
      The input parameter "npiece" describes which term to do.  */
  const BoxArray &ba=SDC.sol[0].boxArray();
  const DistributionMapping &dm=SDC.sol[0].DistributionMap();
  
  const Box& domain_bx = geom.Domain();
  const Real* dx = geom.CellSize();
  ParmParse pp;  
  Real qij;   //  Temp used for quadrature weight
  int numV=0; //  Collects the number of Vcycles
  
  //  Read in the coefficients for A-D-R
  Real a,d,r;
  pp.query("a",a);
  pp.query("d",d);
  pp.query("r",r);
  int Nprob,Lord;
  Real epsilon,k_freq;
  pp.get("Nprob",Nprob);
  pp.get("Lord",Lord);
  pp.get("epsilon",epsilon);
  pp.get("k_freq",k_freq);


  // relative and absolute tolerances for linear solve with MLMG (used to be const)
  Real tol_abs_MG = 1.0e-12;
  pp.query("tol_abs_MG",tol_abs_MG);
  
  Real tol_rel_MG = 1.0e-12;
  pp.query("tol_rel_MG",tol_rel_MG);
  
  // Tolerances for residual in L222(L244)-L4 iteration (used to be const)
  Real tol_abs_res = 1.e-10;    // absolute tolerance on residual
  pp.query("tol_abs_res",tol_abs_res);
  Real tol_rel_res = 1.e-10;    // relative tolerance on residual
  pp.query("tol_rel_res",tol_rel_res);
  
  Real resnorm = 1.e10;    // norm of residual for absolute tolerance
  Real err_rel = 1.e10;  // relative error for relative tolerance
  Real soln_norm = 1.0;
  Real zeroReal = 0.0;

  // Make some space for iteration stuff
  MultiFab corr(ba, dm, 1, 2);
  corr.setVal(1.e10); // so that err_rel in first loop
  MultiFab resid(ba, dm, 1, 0);
  MultiFab eval_storage(ba, dm, 1, 0);
  MultiFab temp_zero(ba, dm, 1, 1);
  temp_zero.setVal(0.0); temp_zero.setBndry(0);
    
  // Set boundary values
  MultiFab bdry_values(ba, dm, 1, 1);
  bdry_values.setVal(0.0); bdry_values.setBndry(0);
  if (Nprob < 3){
    for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
      {          const Box& bx = mfi.validbox();
	fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
			 BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
			 geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&time, &Nprob);
      }
    //  mlabec.fourthOrderBCFill(SDC.sol[sdc_m],bdry_values);
  }
  SDC.sol[sdc_m].FillBoundary(geom.periodicity());
    
  if (npiece == 1)  
    {  // Do diffusion solve
      
        qij = dt*SDC.Qimp[sdc_m-1][sdc_m];  //  Set diffusion scalar in solve
	Real ascalar = 1.0;
        mlabec.setScalars(ascalar, d*qij);  //  Set "mass matrix"
        
        int res_iter=0;   //  Initialize residual loop counter

	int max_res_iter = 50;  //  Set max residual loops
        pp.query("max_res_iter",max_res_iter);

	int max_Vcycle=3;  // Set max number of Vcycles per residual solve
	pp.query("max_Vcycle",max_Vcycle);
        
        while ((resnorm > tol_abs_res) & (res_iter <=max_res_iter))  //  Loop over residual solves	  
	  {
	    //  First compute the residual by calling the diffusion operator and subtracting rhs
	    //  Compute the diffusion operator
            for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
            {
	      const Box& bx = mfi.validbox();
	      SDC_feval_F(BL_TO_FORTRAN_BOX(bx),
			  BL_TO_FORTRAN_BOX(domain_bx),
			  BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
			  BL_TO_FORTRAN_ANYD(flux[0][mfi]),
			  BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)
			  BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
			  BL_TO_FORTRAN_ANYD(eval_storage[mfi]),
			  dx,&a,&d,&r,
			  BL_TO_FORTRAN_ANYD(face_bcoef[0][mfi]),
			  BL_TO_FORTRAN_ANYD(face_bcoef[1][mfi]),
			  BL_TO_FORTRAN_ANYD(prod_stor[0][mfi]),
			  BL_TO_FORTRAN_ANYD(prod_stor[1][mfi]),
			  &npiece, &time, &epsilon, &k_freq, &Nprob, &Lord);			  

	      
            }
            //Rescale eval_storage to make resid.
            resid.setVal(0.0);
	    MultiFab::Copy(resid, rhs,0, 0, 1, 0);	    
	    MultiFab::Saxpy(resid,qij,eval_storage,0,0,1,0);
	    MultiFab::Saxpy(resid,-1.0,SDC.sol[sdc_m],0,0,1,0);

            //corrnorm=eval_storage.norm0();
            //amrex::Print() << "iter " << res_iter << ",  Diffusion operator norm " << corrnorm << "\n";
            
	    resnorm = resid.norm0();
	    err_rel = corr.norm0()/soln_norm;
	    
	    if(resnorm <= tol_abs_res ){
	      amrex::Print() << "---+--- Reached absolute residual tolerance in L2-L4: res_iter=" << res_iter << ",  residual norm=" << resnorm << ", relative error=" << err_rel << "\n";
	      break;
	    }
	    if(err_rel <= tol_rel_res ){
              amrex::Print() << "---+--- Reached relative resdiual tolerance in L2-L4: res_iter=" << res_iter << ",  residual norm=" << resnorm << ", relative error=" << err_rel << "\n";
              break;
	    }

	    //  Set some BC (since this is a correction, they are zero
	    mlabec.setLevelBC(0,&temp_zero); 
            corr.setVal(0.0);

	    //  Do the multigrid solve for residual
	    mlmg.setFixedIter(max_Vcycle);                
	    mlmg.setVerbose(0);                
	    mlabec.prepareForSolve();
	    mlmg.solve({&corr}, {&resid}, tol_rel_MG, tol_abs_MG);
	    numV=mlmg.getNumIters();
	    tot_Vcycle+=numV;
	    tot_res_iter++;
	    amrex::Print() << "---+--- Residual iter=" << res_iter << ",  residual norm=" << resnorm << ",  numV=" << numV << "\n";	    
	    soln_norm = SDC.sol[sdc_m].norm0();
	    
	    //  Add correction to solution
	    MultiFab::Saxpy(SDC.sol[sdc_m],1.0,corr,0,0,1,0);
	    
	    // Set boundaries on new solution
            if (Nprob<3){mlabec.fourthOrderBCFill(SDC.sol[sdc_m],bdry_values);}
            SDC.sol[sdc_m].FillBoundary(geom.periodicity());
            ++res_iter;
	  }
    }
  else
    {  // Do reaction solve  y - qij*y*(1-y)*(y-1/2) = rhs

      //  make a flag to change how the reaction is done
      int nflag=1;  // Lazy approximation
      
      qij = r*dt*SDC.Qimp[sdc_m-1][sdc_m];	            
      for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
	{
	  const Box& bx = mfi.validbox();
	  SDC_fcomp_reaction_F(BL_TO_FORTRAN_BOX(bx),
			       BL_TO_FORTRAN_BOX(domain_bx),
			       BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
			       BL_TO_FORTRAN_ANYD(rhs[mfi]),		      
			       BL_TO_FORTRAN_ANYD(SDC.f[2][sdc_m][mfi]),
			       &qij,&nflag); 
	}
    }

}




