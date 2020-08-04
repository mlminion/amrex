#include "AMReX_SDCstruct.H"


SDCstruct::SDCstruct(int Nnodes_in,int Npieces_in, MultiFab& sol_in)
{
  
  Nnodes=Nnodes_in;
  Npieces=Npieces_in;       
  Ncomp=sol_in.nComp();
  
  qnodes= new Real[Nnodes];
  Qall= new Real[4*(Nnodes-1)*Nnodes];  
  Nflags= new int[Nnodes];

  Qgauss.resize(Nnodes-1, Vector<Real>(Nnodes));
  Qexp.resize(Nnodes-1, Vector<Real>(Nnodes));
  Qimp.resize(Nnodes-1, Vector<Real>(Nnodes));
  QLU.resize(Nnodes-1, Vector<Real>(Nnodes));  

  //  Make the quadrature tables
  SDC_quadrature(&qtype, &Nnodes, &Nnodes,qnodes,Nflags, &Qall[0]);  

  //  Load the quadrature nodes into their spots
  for ( int j = 0; j < Nnodes-1; ++j)
    for ( int k = 0; k < Nnodes; ++k)
      {
	Qgauss[j][k]=  Qall[0*(Nnodes-1)*(Nnodes) +j*(Nnodes) + k ];
	Qexp[j][k]=    Qall[1*(Nnodes-1)*(Nnodes) +j*(Nnodes) + k ];
	Qimp[j][k]=    Qall[2*(Nnodes-1)*(Nnodes) +j*(Nnodes) + k ];
	QLU[j][k]=     Qall[3*(Nnodes-1)*(Nnodes) +j*(Nnodes) + k ];			
      }

  //  Resize the storage
  sol.resize(Nnodes);
  res.resize(Nnodes);
  f.resize(Npieces);
  if (Npieces == 3) 
    Ithree.resize(Nnodes);

  //  Assign  geomety and multifab info
  const BoxArray &ba=sol_in.boxArray();
  const DistributionMapping &dm=sol_in.DistributionMap();
  const int Nghost=sol_in.nGrow();

  for (auto& v : f) v.resize(Nnodes);  
  for (int sdc_m = 0; sdc_m < Nnodes; sdc_m++)
    {
      sol[sdc_m].define(ba, dm, Ncomp, Nghost);
      res[sdc_m].define(ba, dm, Ncomp, 0);
      for (int i = 0; i < Npieces; i++)
	{
	  f[i][sdc_m].define(ba, dm, Ncomp, 0);
	}
      if (Npieces == 3)
	Ithree[sdc_m].define(ba, dm, Ncomp, 0);      
    }
  
}

void SDCstruct::SDC_integrals(Real dt)
{

  Real qij;
  Real qijE,qijI;

  // Compute the quadrature terms from last iteration
  for (int sdc_m = 0; sdc_m < Nnodes-1; sdc_m++)
    {
      //    MultiFab:Setval(res[sdc_m],0.0);
      res[sdc_m]=0.0;      
      if (Npieces == 3)
	Ithree[sdc_m]=0.0;

      for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
	{
	  qij = dt*(Qgauss[sdc_m][sdc_n]);	      
	  MultiFab::Saxpy(res[sdc_m],qij,f[0][sdc_n],0, 0, 1, 0);
	  qij = dt*(Qgauss[sdc_m][sdc_n]);	      	  
	  MultiFab::Saxpy(res[sdc_m],qij,f[1][sdc_n],0, 0, 1, 0);
	  if (Npieces == 3)
	    {
	      qij = dt*(Qgauss[sdc_m][sdc_n]);  // leave off -dt*Qtil and add it later
	      MultiFab::Saxpy(res[sdc_m],qij,f[2][sdc_n],0, 0, 1, 0);	      
	      //	      res[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
	    }
	  
	}
      /*
      for ( MFIter mfi(res[sdc_m]); mfi.isValid(); ++mfi )
	{
	  //	  res[sdc_m][mfi].setVal(0.0);
	  if (Npieces == 3)
	    Ithree[sdc_m][mfi].setVal(0.0);
	  for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
	    {
	      qij = dt*(Qgauss[sdc_m][sdc_n]-Qexp[sdc_m][sdc_n]);	      
	      //	      res[sdc_m][mfi].saxpy(qij,f[0][sdc_n][mfi]);  // Explicit part

	      qij = dt*(Qgauss[sdc_m][sdc_n]-Qimp[sdc_m][sdc_n]);	      
	      //	      res[sdc_m][mfi].saxpy(qij,f[1][sdc_n][mfi]);  // Implicit part
	    }
	  if (Npieces == 3)
	    {
	      for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
		{ //  // MISDC pieces
		  qij = dt*(Qgauss[sdc_m][sdc_n]);  // leave off -dt*Qtil and add it later		  
		  res[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
		  // Compute seperate integral for f_3 piece		  
		  qij = -dt*(Qimp[sdc_m][sdc_n]);  
		  Ithree[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
		}
	    }
	}
      */
    }
}
Real SDCstruct::SDC_check_resid(MultiFab& res_m, Real dt,int m)
{
  //  Compute the residual at node m and return the norm
  //  When this is called, the integrals must have been previously set
  //  by a call to SDC_integrals
  Real norm;
  
  MultiFab::Copy(res_m,sol[0], 0, 0, 1, 0);
  MultiFab::Add(res_m,res[m-1], 0, 0, 1, 0);
  MultiFab::Subtract(res_m,sol[m], 0, 0, 1, 0);

  norm = res_m.norm0()/dt;  
  return norm;
  
}
void SDCstruct::SDC_integrals_diff(Real dt)
{

  Real qij;
  Real qijE,qijI;

  // Compute the quadrature terms from last iteration
  for (int sdc_m = 0; sdc_m < Nnodes-1; sdc_m++)
    {
      for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
	{
	  qij = dt*(-Qexp[sdc_m][sdc_n]);	      
	  MultiFab::Saxpy(res[sdc_m],qij,f[0][sdc_n],0, 0, 1, 0);
	  qij = dt*(-Qimp[sdc_m][sdc_n]);	      	  
	  MultiFab::Saxpy(res[sdc_m],qij,f[1][sdc_n],0, 0, 1, 0);
	  if (Npieces == 3)
	    {
	      // Compute seperate integral for f_3 piece		  
	      qij = -dt*(Qimp[sdc_m][sdc_n]);  
	      //	      Ithree[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
	      MultiFab::Saxpy(Ithree[sdc_m],qij,f[2][sdc_n],0, 0, 1, 0);	      
	    }
	  
	}
      /*
      for ( MFIter mfi(res[sdc_m]); mfi.isValid(); ++mfi )
	{
	  //	  res[sdc_m][mfi].setVal(0.0);
	  if (Npieces == 3)
	    Ithree[sdc_m][mfi].setVal(0.0);
	  for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
	    {
	      qij = dt*(Qgauss[sdc_m][sdc_n]-Qexp[sdc_m][sdc_n]);	      
	      //	      res[sdc_m][mfi].saxpy(qij,f[0][sdc_n][mfi]);  // Explicit part

	      qij = dt*(Qgauss[sdc_m][sdc_n]-Qimp[sdc_m][sdc_n]);	      
	      //	      res[sdc_m][mfi].saxpy(qij,f[1][sdc_n][mfi]);  // Implicit part
	    }
	  if (Npieces == 3)
	    {
	      for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
		{ //  // MISDC pieces
		  qij = dt*(Qgauss[sdc_m][sdc_n]);  // leave off -dt*Qtil and add it later		  
		  res[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
		  // Compute seperate integral for f_3 piece		  
		  qij = -dt*(Qimp[sdc_m][sdc_n]);  
		  Ithree[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
		}
	    }
	}
      */
    }
}
void SDCstruct::SDC_rhs_integrals(Real dt)
{

  Real qij;
  Real qijE,qijI;

  // Compute the quadrature terms from last iteration
  for (int sdc_m = 0; sdc_m < Nnodes-1; sdc_m++)
    {
      //    MultiFab:Setval(res[sdc_m],0.0);
      res[sdc_m]=0.0;      
      if (Npieces == 3)
	Ithree[sdc_m]=0.0;

      for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
	{
	  qij = dt*(Qgauss[sdc_m][sdc_n]-Qexp[sdc_m][sdc_n]);	      
	  MultiFab::Saxpy(res[sdc_m],qij,f[0][sdc_n],0, 0, 1, 0);
	  qij = dt*(Qgauss[sdc_m][sdc_n]-Qimp[sdc_m][sdc_n]);	      	  
	  MultiFab::Saxpy(res[sdc_m],qij,f[1][sdc_n],0, 0, 1, 0);
	  if (Npieces == 3)
	    {
	      qij = dt*(Qgauss[sdc_m][sdc_n]);  // leave off -dt*Qtil and add it later
	      MultiFab::Saxpy(res[sdc_m],qij,f[2][sdc_n],0, 0, 1, 0);	      
	      //	      res[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
	      // Compute seperate integral for f_3 piece		  
	      qij = -dt*(Qimp[sdc_m][sdc_n]);  
	      //	      Ithree[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
	      MultiFab::Saxpy(Ithree[sdc_m],qij,f[2][sdc_n],0, 0, 1, 0);	      
	    }
	  
	}
      /*
      for ( MFIter mfi(res[sdc_m]); mfi.isValid(); ++mfi )
	{
	  //	  res[sdc_m][mfi].setVal(0.0);
	  if (Npieces == 3)
	    Ithree[sdc_m][mfi].setVal(0.0);
	  for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
	    {
	      qij = dt*(Qgauss[sdc_m][sdc_n]-Qexp[sdc_m][sdc_n]);	      
	      //	      res[sdc_m][mfi].saxpy(qij,f[0][sdc_n][mfi]);  // Explicit part

	      qij = dt*(Qgauss[sdc_m][sdc_n]-Qimp[sdc_m][sdc_n]);	      
	      //	      res[sdc_m][mfi].saxpy(qij,f[1][sdc_n][mfi]);  // Implicit part
	    }
	  if (Npieces == 3)
	    {
	      for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
		{ //  // MISDC pieces
		  qij = dt*(Qgauss[sdc_m][sdc_n]);  // leave off -dt*Qtil and add it later		  
		  res[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
		  // Compute seperate integral for f_3 piece		  
		  qij = -dt*(Qimp[sdc_m][sdc_n]);  
		  Ithree[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
		}
	    }
	}
      */
    }
}

void SDCstruct::SDC_rhs_k_plus_one(MultiFab& sol_new, Real dt,int sdc_m)
{
  //  Compute the rhs terms for the implicit solve
  Real qij;
  
  //  Copy first the initial value
  MultiFab::Copy(sol_new,sol[0], 0, 0, 1, 0);
  MultiFab::Add(sol_new,res[sdc_m], 0, 0, 1, 0);

  for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
    {
      	  qij = dt*Qexp[sdc_m][sdc_n];
	  MultiFab::Saxpy(sol_new,qij,f[0][sdc_n],0, 0, 1, 0);
      	  qij = dt*Qimp[sdc_m][sdc_n];
	  MultiFab::Saxpy(sol_new,qij,f[1][sdc_n],0, 0, 1, 0);
    }
  //  for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
  //    {
  //      	  qij = dt*Qexp[sdc_m][sdc_n];
  //	  MultiFab::Saxpy(sol_new,qij,f[0][sdc_n],0, 0, 1, 0);
  //    }
  
  //	{
  //
  //	  sol_new[mfi].saxpy(qij,f[0][sdc_n][mfi]);
  //	}
  //  for ( MFIter mfi(sol_new); mfi.isValid(); ++mfi )
  //    {
      //      sol_new[mfi].Add(res[sdc_m][mfi]);
  //      for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
  //	{
  //	  qij = dt*Qexp[sdc_m][sdc_n];
  //	  sol_new[mfi].saxpy(qij,f[0][sdc_n][mfi]);
  //	}
  //      for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
  //	{
  //	  qij = dt*Qimp[sdc_m][sdc_n];
  //	  sol_new[mfi].saxpy(qij,f[1][sdc_n][mfi]);
  //	}
  //  }
  
}
void SDCstruct::SDC_rhs_misdc(MultiFab& sol_new, Real dt,int sdc_m)
{
  //  Add the terms to the rhs before the second implicit solve
  Real qij;
  
  for ( MFIter mfi(sol_new); mfi.isValid(); ++mfi )
    {
      sol_new[mfi].saxpy(1.0,Ithree[sdc_m][mfi]);
      for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
	{
	  qij = dt*Qimp[sdc_m][sdc_n];
	  sol_new[mfi].saxpy(qij,f[2][sdc_n][mfi]);
	}
    }
  
}
