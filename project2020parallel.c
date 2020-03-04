/**********************************************************************
 *
 * Diaforiki exiswsi - Provlima arxikwn timwn - Parallel
 * Processor 2.6GHz Intel Core i5
 *
 *********************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>

int main(){
  int i,j,k,Lx=10,Ly=10,Nx=40,Ny=40,T1=20,c=1,nt;
  double h,dx,dy,dt,cx,cy,X[40],Y[40],I[40][40],V[40][40],Un[40][40],Uan,Uan5,Uan10,Uan20,Unplus1[40][40],Uanal[40][40],
    fTimeStart,fTimeEnd;

  // Record start time
  fTimeStart = omp_get_wtime();

  h=(double)Lx/Nx;
  dt=h*h/(2.*c);
  dx=(double)Lx/(Nx-1);
  dy=(double)Ly/(Ny-1);
  nt=(double)T1/dt;
  double t[nt];
  cx=(double)c*dt/dx; /*xwris tetragwno*/
  cy=(double)c*dt/dy;

  FILE *data0;
  data0=fopen("data0.txt","w");
  FILE *data5;
  data5=fopen("data5.txt","w");
  FILE *data10;
  data10=fopen("data10.txt","w");
  FILE *data20;
  data20=fopen("data20.txt","w");
  FILE *data0anal;
  data0anal=fopen("data0anal.txt","w");
  FILE *data5anal;
  data5anal=fopen("data5anal.txt","w");
  FILE *data10anal;
  data10anal=fopen("data10anal.txt","w");
  FILE *data20anal;
  data20anal=fopen("data20anal.txt","w");
  
  X[0]=0;
  Y[0]=0;
  omp_set_dynamic(0);
#pragma omp parallel private(i,j,k) shared(h,nt,t,dt,Nx,Ny,Lx,Ly,I,V,X,Y,dx,dy,Uan,Uan5,Uan10,Uan20,Uanal,Un,Unplus1,c,cx,cy,data0,data0anal,data5,data5anal,data10,data10anal,data20,data20anal) default(none)
  {    
#pragma omp for
    for(i=1;i<40;i++){
      X[i]=X[0]+i*h;
      Y[i]=Y[0]+i*h;
    }

#pragma omp for
    for (k=0; k<=nt; k++){t[k]=k*dt;} //Orismos xronou

#pragma omp for 
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	I[i][j]=X[i]*(Lx-X[i])*Y[j]*(Ly-Y[j]);
	V[i][j]=X[i]*(Lx-X[i])*Y[j]*(Ly-Y[j])/2.;
	k=0;
	Uanal[i][j]=X[i]*(Lx-X[i])*Y[j]*(Ly-Y[j])*(1.+t[k]/2.);
      }
    }

#pragma omp single
    {
      Uan=X[Nx/2]*(Lx-X[Nx/2])*Y[Ny/2]*(Ly-Y[Ny/2])*(1+0./2);
      printf("Time =%f\t, Calculated value at (Nx/2,Ny/2) = %f,  Analytical value = %f\n",
	     k*dt,I[Nx/2][Ny/2],Uan);	
    }
    
#pragma omp for
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	fprintf(data0,"%f ",I[i][j]);
	fprintf(data0anal,"%f ",Uanal[i][j]);
      }
      fprintf(data0,"\n");  // thn t=0 exoume mono I(x,y)
      fprintf(data0anal,"\n");
    }
    
    //Finite Difference
#pragma omp for 
    for(j=1;j<Ny-1;j++){
      for(i=1;i<Nx-1;i++){
	Un[i][j]=I[i][j]+dt*V[i][j]+(I[i+1][j]-2.*I[i][j]+I[i-1][j])*cx*cx/2.+ \
	  (I[i][j+1]-2.*I[i][j]+I[i][j-1])*cy*cy/2.+			\
	  dt*dt*c*c*(1.+1.0/2.)*(Y[j]*(Ly-Y[j])+X[i]*(Lx-X[i]));
	// thn t=1 exoume thn Un
      }
    }
    
//xronos apo t=2 ews Nt,opou Nt=640,t1=5 simenei 5/dt=160,10/dt=321,20/dt=640
    for(k=2;k<=nt;k++){
#pragma omp for schedule(static)
      for(j=1;j<Ny-1;j++){
	for(i=1;i<Nx-1;i++){
	  Unplus1[i][j]=-I[i][j]+2*Un[i][j]+(Un[i+1][j]-2*Un[i][j]+Un[i-1][j])*cx*cx+ \
	    (Un[i][j+1]-2*Un[i][j]+Un[i][j-1])*cy*cy+			\
	    dt*dt*2*c*c*(1.+t[k]/2.)*(Y[j]*(Ly-Y[j])+X[i]*(Lx-X[i]));
	  Uanal[i][j]=X[i]*(Lx-X[i])*Y[j]*(Ly-Y[j])*(1.+t[k]/2.);}}
      
      
      if (k==5){
#pragma omp single
	{
	  Uan5=X[Nx/2]*(Lx-X[Nx/2])*Y[Ny/2]*(Ly-Y[Ny/2])*(1+t[k]/2.);
	  printf("Time =%f\t, Calculated value at (Nx/2,Ny/2) = %f,  Analytical value = %f\n",
		 k*dt,I[Nx/2][Ny/2],Uan5);	
	}	
      }
      if (k==10){
#pragma omp single
	{
	  Uan10=X[Nx/2]*(Lx-X[Nx/2])*Y[Ny/2]*(Ly-Y[Ny/2])*(1+t[k]/2.);
	  printf("Time =%f\t, Calculated value at (Nx/2,Ny/2) = %f,  Analytical value = %f\n",
		 k*dt,I[Nx/2][Ny/2],Uan10);	
	}	
      }

      if (k==20){
#pragma omp single
	{
	  Uan20=X[Nx/2]*(Lx-X[Nx/2])*Y[Ny/2]*(Ly-Y[Ny/2])*(1+t[k]/2.);
	  printf("Time =%f\t, Calculated value at (Nx/2,Ny/2) = %f,  Analytical value = %f\n",
		 k*dt,I[Nx/2][Ny/2],Uan20);
	}
      }
      
//swap arrays gia na leitourgisei o algorithmos alliws thelei 3sdiastato pinaka
#pragma omp for
      for(j=1;j<Ny-1;j++){
      for(i=1;i<Nx-1;i++){
	I[i][j]=Un[i][j];
      }
      }
      
#pragma omp for
      for(j=1;j<Ny-1;j++){
	for(i=1;i<Nx-1;i++){
	  Un[i][j]=Unplus1[i][j];
	}
      }

#pragma omp for
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
#pragma omp critical
	  {
	    if(k==160){
	      fprintf(data5,"%f ",Unplus1[i][j]);
	      fprintf(data5anal,"%f ",Uanal[i][j]);}
	    else if(k==321){
	      fprintf(data10,"%f ",Unplus1[i][j]);
	      fprintf(data10anal,"%f ",Uanal[i][j]);}
	    else if(k==640){
	      fprintf(data20,"%f ",Unplus1[i][j]);
	      fprintf(data20anal,"%f ",Uanal[i][j]);}
	  }
	}
	if(k==160){
	  fprintf(data5,"\n ");
	  fprintf(data5anal,"\n ");}
	else if(k==321){
	  fprintf(data10,"\n ");
	  fprintf(data10anal,"\n ");}
	else if(k==640){
	  fprintf(data20,"\n ");
	  fprintf(data20anal,"\n ");}
      }
      
    } //for xronou

  }//parallel
  
  fclose(data0);
  fclose(data0anal);
  fclose(data5);
  fclose(data5anal);
  fclose(data10);
  fclose(data10anal);
  fclose(data20);
  fclose(data20anal);
  
  // Record end time
  fTimeEnd = omp_get_wtime();

// Print elapsed time
  printf("wall clock time= %.20f\n", fTimeEnd - fTimeStart);

  return 0; 
}//main
