//
// Stochastic simulacion of an exotic tree invasion
//

#include "bgi.hpp"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
extern "C"
{
#include "\cpp\randlib\src\randlib.h"
}


#define M_PI		3.14159265358979323846

#include "smattpl.h"
#include <RWFile.h>


struct parameters
{
    char gr;
    long rndSeed;
    int DimX;
    int DimY;
    int tFinal;
    bool fullStoch;
    int writeMap;
    bool pomac;
    double g;
    double dj0;
    double dj1;
    double alfa;
    double b;
    double da;
    double densIni;
    char fName[20];
};

int readParms(char * pName,parameters &p);

double StochasticInvasion(parameters p,simplmat <long> &nTotal,simplmat <long> &NTotal,
				simplmat <long> &nOcupado,simplmat <long> &NOcupado);

double DeterministicInvasion(parameters p,simplmat <long> &nTotal,simplmat <long> &NTotal,
				simplmat <long> &nOcupado,simplmat <long> &NOcupado);
				
int MuestreoPorCaminataAlAzar(simplmat <double> &NN,simplmat <double> &nn,
				simplmat <double> &prome);
				
int ReadPomacParms(ifstream &inLineFile,parameters &p );				
				
int main(int argc, char * argv[])
{
	int i;
	parameters ap;

	long seed1=0,seed2=0; 

	simplmat <long> nTotal;
	simplmat <long> NTotal;
	simplmat <long> nOcupado;
	simplmat <long> NOcupado;
	bool priPro=false;
	
   	if( argc==2 )
   	{
        if( readParms(argv[1],ap) )
        {
        	cerr << "bad parms file " << argv[1] << endl;
        	exit(1);
        }
		seed1 = ap.rndSeed;
		if(seed1==-1)
			seed1=time(0);
		seed2 = (seed1+1);

		setall(seed1,seed2);

		nTotal.resize(ap.tFinal+1);
		nTotal.fill(0);
		NTotal.resize(ap.tFinal+1);
		NTotal.fill(0);
		NOcupado.resize(ap.tFinal+1);
		NOcupado.fill(0);
		nOcupado.resize(ap.tFinal+1);
		nOcupado.fill(0);
		double timeOcupa;
		
		if( ap.pomac )
		{
			ifstream pomacIn;
			while( ReadPomacParms(pomacIn,ap) )
			    if( ap.fullStoch )
				    timeOcupa = StochasticInvasion(ap,nTotal,NTotal,nOcupado,NOcupado);
			    else
				    timeOcupa = DeterministicInvasion(ap,nTotal,NTotal,nOcupado,NOcupado);
		}
		else
		{
			if( ap.fullStoch )
				timeOcupa = StochasticInvasion(ap,nTotal,NTotal,nOcupado,NOcupado);
			else
				timeOcupa = DeterministicInvasion(ap,nTotal,NTotal,nOcupado,NOcupado);
		}	
		
		if(!ap.pomac)
		{
			ostringstream fname;
			fname << ap.fName << ".out"<< ends;

            
			fstream datafile;
			datafile.open(fname.str().c_str(),ios::in);
			if(!datafile )
				priPro=true;
			datafile.close();
			datafile.clear();
			datafile.open(fname.str().c_str(), ios::out | ios::app);
			if( !datafile )
				{
					cerr << "Cannot open file: " << ap.fName << endl;
					return 0;
				}
				if( priPro )
				{
					datafile << "DimX\tDimY\ttFinal\tg\tdj0\tdj1\talfa\tb\tda\tdensIni\tTime\tJuv\tAd\tJuvOcupa\tAdOcupa" << endl;				
					priPro = false;
				}

			for(i=0; i<ap.tFinal+1; i++)
			{
//				if(nOcupado(i)==0 && NOcupado(i)==0)
//					break;
				datafile << ap.DimX << "\t"
					<< ap.DimY << "\t"
					<< ap.tFinal << "\t"
					<< ap.g << "\t"
					<< ap.dj0 << "\t"
					<< ap.dj1 << "\t"
					<< ap.alfa << "\t"
					<< ap.b << "\t"
					<< ap.da << "\t"
					<< ap.densIni << "\t"
					<< i << "\t"
					<< nTotal(i) << "\t"
					<< NTotal(i) << "\t"
					<< nOcupado(i) << "\t"
					<< NOcupado(i) << endl;
			}	

			datafile << ap.DimX << "\t"
					<< ap.DimY << "\t"
					<< ap.tFinal << "\t"
					<< ap.g << "\t"
					<< ap.dj0 << "\t"
					<< ap.dj1 << "\t"
					<< ap.alfa << "\t"
					<< ap.b << "\t"
					<< ap.da << "\t"
					<< ap.densIni << "\t"
					<< "Time50%\t" 
					<< timeOcupa << endl;

		}

		return(0);

	}
	else
	{
        cerr << "Syntaxis: inva parms.dat" << endl;
    	exit(1);
    }
        
	return 0;
}


int readParms(char * pName,parameters &p)
	{
   	char buff[50];
	ifstream parms(pName);
   	if( !parms )
		{
      	cerr << "Can't open parms file" << endl;
      	return(1);
      	}

	parms >> buff;
	if(strcmp(buff,"Gr")==0)
    	parms >> p.gr;
    else
    	return(1);
	p.gr = toupper( p.gr );

	parms >> buff;
    if(strcmp(buff,"RndSeed")==0)
    	parms >> p.rndSeed;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"DimX")==0)
    	parms >> p.DimX;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"DimY")==0)
    	parms >> p.DimY;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"tFinal")==0)
    	parms >> p.tFinal;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"fullStoch")==0)
    	parms >> p.fullStoch;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"writeMap")==0)
    	parms >> p.writeMap;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"pomac")==0)
    	parms >> p.pomac;
    else
    	return(1);
    	
	parms >> buff;
    if(strcmp(buff,"g")==0)
    	parms >> p.g;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"dj0")==0)
    	parms >> p.dj0;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"dj1")==0)
    	parms >> p.dj1;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"alfa")==0)
    	parms >> p.alfa;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"b")==0)
    	parms >> p.b;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"da")==0)
    	parms >> p.da;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"densIni")==0)
    	parms >> p.densIni;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"fName")==0)
    	parms >> p.fName;
    else
    	return(1);

    parms.close();

    return(0);
    }



double StochasticInvasion(parameters p,simplmat <long> &nTotal,simplmat <long> &NTotal,
				simplmat <long> &nOcupado,simplmat <long> &NOcupado)
{
	int privez=1;
	
	int i,j,priVez=1;
    double rnd;

	simplmat <long> N(p.DimX,p.DimY);
	simplmat <long> n(p.DimX,p.DimY);
//	simplmat <long> NAnt1(p.Dim,p.Dim);
	
	N.fill(0);
	n.fill(0);
	
#ifdef GRAPHICS
	if(p.gr=='S')
		{
		IGraph(p.DimX,p.DimY);
		ReadIdrisiPalette("Idris256.smp");
		}
#endif

	//
	// Condiciones iniciales fila de 4 celdas con 10 adultos sobre último X
	//
    //N(p.DimX-1,p.DimY/2)=p.densIni;
    //N(p.DimX-1,(p.DimY/2)-1)=p.densIni;
    //N(p.DimX-1,(p.DimY/2)-2)=p.densIni;
    //N(p.DimX-1,(p.DimY/2)+1)=p.densIni;

	//
	// Condiciones iniciales (10 adultos en las cuatro celdas centrales)
	//
    N(p.DimX/2,p.DimY/2)=p.densIni;
    N(p.DimX/2,(p.DimY/2)-1)=p.densIni;
    N((p.DimY/2)-1,p.DimY/2)=p.densIni;
    N((p.DimY/2)-1,(p.DimY/2)-1)=p.densIni;
    
	//
	// Condiciones iniciales 50% al azar 
	//
	//	for(i=0; i<p.Dim; i++)
	//		for(j=0; j<p.Dim; j++)
	//    	{
	//     		rnd = ranf();
	//			n(i,j) = (rnd > 0.5) ? 1 : 0;
	//     	}

   	long totalJuv = 0,totalAd = 0;

    double time=0,totalRate=0,eventRate,time50=0;
    int dTime=0,dTimeAnt=0;
    int jPlus1,iPlus1,iMinus1,jMinus1;
    long juvOcupado=0,adOcupado=0;

	int salida=0;

	while( 1 ) 
	{
        totalJuv=0;
        totalAd=0;
        totalRate=0;
        juvOcupado=0;
        adOcupado=0;

      	for( i=0; i<p.DimX; i++)
      		for( j=0; j<p.DimY; j++)
			{
				totalJuv += n(i,j);
				totalAd += N(i,j);

				if( N(i,j)>0 )
					adOcupado++;

				if( n(i,j)>0 )
					juvOcupado++;


				totalRate += ( p.g + p.dj0 + p.dj1*N(i,j) ) * n(i,j) + ((1-p.alfa)* p.b + p.da) * N(i,j) ;

#ifdef NOWRAP
                iPlus1  = i + 1;
                jPlus1  = j + 1;
                iMinus1 = i - 1;
                jMinus1 = j - 1;

#else
                iPlus1  = (p.DimX + i + 1) % p.DimX;
                jPlus1  = (p.DimY + j + 1) % p.DimY;
                iMinus1 = (p.DimX + i - 1) % p.DimX;
                jMinus1 = (p.DimY + j - 1) % p.DimY;
#endif

				totalRate += p.alfa/8 * p.b * ( 									
									N(iPlus1,jPlus1)+N(iPlus1,j)+N(iPlus1,jMinus1)+
									N(i,jPlus1)+N(i,jMinus1)+
									N(iMinus1,jPlus1)+N(iMinus1,j)+N(iMinus1,jMinus1)
									);
           	}
		//
		// Hacer salida 
		//
		if( dTime>dTimeAnt || priVez || salida)
		{
			nTotal(dTime)= totalJuv;
			NTotal(dTime)= totalAd;
			nOcupado(dTime) = juvOcupado;
			NOcupado(dTime) = adOcupado;
			dTimeAnt = dTime;

			if(adOcupado >= p.DimX*p.DimY/2 && time50==0)
			{
				time50= time;       // Guarda el tiempo al llegar al 50% de ocupación
				salida=1;           //  Termina cuando llegaba al 50% de la superficie ocupada
			}

			// En forma de mapa
			//

            if( p.writeMap>0 && p.writeMap == dTime )
            {
            	RWFile file;
				ostringstream fname;
				fname << p.fName << dTime << ".map"<< ends;
			
				if(!file.WriteSeed(fname.str().c_str(),N, "AD"))
                    cerr << "No pudo escribir salida AD" << endl;
				if(!file.WriteSeed(fname.str().c_str(),n, "JU"));
                    cerr << "No pudo escribir salida JU" << endl;
			}
			
			cerr << dTime<< "\t" << totalJuv << "\t" << totalAd	<< "\t" << adOcupado << endl;
			
            if( p.pomac && salida)
            {    
                
                ostringstream pname;
                bool priPro=false;
                
                pname << p.fName << "Pomac.out" << ends;
                
	        	fstream pomacfile;
				pomacfile.open(pname.str().c_str(),ios::in);
				if(!pomacfile )
					priPro=true;
				pomacfile.close();
				pomacfile.clear();

        		pomacfile.open(pname.str().c_str(), ios::out | ios::app);
		        if( !pomacfile )
    			{
				    cerr << "Cannot open file: " << pname << endl;
				    return 0;
	    		}
				if( priPro )
				{
					pomacfile << "g\tdj0\tdj1\talfa\tb\tda\tdensIni\tadOcupado\ttime50\tVel" << endl;
					priPro = false;
				}

                pomacfile << p.g << "\t"
    			        << p.dj0 << "\t"
            			<< p.dj1 << "\t"
            			<< p.alfa << "\t"
        				<< p.b << "\t"
        				<< p.da << "\t"
        				<< p.densIni << "\t"
						<< adOcupado << "\t"
                        << time50 << "\t"
                        << sqrt(double(adOcupado))/time50*100 << endl;
			}   

          	priVez=0;
			/*
	      	for( j=0; j<p.DimY; j++)
			{
    	  		for( i=0; i<p.DimX; i++)
					cout << n(i,j) << "\t";
				cout << -1 << "\t" << -1 << "\t";
    	  		for( i=0; i<p.DimX; i++)
					cout << N(i,j) << "\t";
				cout << endl;
			}
   	  		for( i=0; i<p.DimX; i++)
				cout << -1 << "\t" << -1 << "\t";
			cout << -1 << "\t" << -1 << endl;
			*/
#ifdef DELAYREGE
			//
			// Guardar valores anteriores
			//
	      	for( i=0; i<p.Dim; i++)
    	  		for( j=0; j<p.Dim; j++)
					NAnt1(i,j)=N(i,j);
#endif

#ifdef GRAPHICS
			if(p.gr=='S')
	      		for( j=0; j<p.DimY; j++)
    	  			for( i=0; i<p.DimX; i++)
         				PPix(i,j,N(i,j)%256);
#endif
			if(salida)
				return time;

			if( dTime==p.tFinal )
				return time50;
		}


        rnd =ranf();
        eventRate = 0;
		if( totalRate==0)
          return 0;
			
      	for( i=0; i<p.DimX; i++)
      		for( j=0; j<p.DimY; j++)
			{
				//
             	// Pasa de Juvenil a Adulto
				//
				eventRate +=  p.g*n(i,j)/totalRate;
                if(eventRate>=rnd)
                	{
                  	n(i,j)--;
                    N(i,j)++;
                    goto finEval;
                    }
				else

				//
				// Mortalidad de Juveniles
				//
				eventRate +=  (p.dj0+p.dj1*N(i,j))*n(i,j)/totalRate;
                if(eventRate>=rnd)
                	{
                  	n(i,j)--;
                    goto finEval;
                    }


				//
				// Regeneracion
				//
#ifdef NOWRAP
                iPlus1  = i + 1;
                jPlus1  = j + 1;
                iMinus1 = i - 1;
                jMinus1 = j - 1;

#else
                iPlus1  = (p.DimX + i + 1) % p.DimX;
                jPlus1  = (p.DimY + j + 1) % p.DimY;
                iMinus1 = (p.DimX + i - 1) % p.DimX;
                jMinus1 = (p.DimY + j - 1) % p.DimY;
#endif

				eventRate +=  ((1-p.alfa)*p.b*N(i,j)+  
								p.alfa/8*p.b*(
									N(iPlus1,jPlus1)+N(iPlus1,j)+N(iPlus1,jMinus1)+
									N(i,jPlus1)+N(i,jMinus1)+
									N(iMinus1,jPlus1)+N(iMinus1,j)+N(iMinus1,jMinus1)))/totalRate;
                if(eventRate>=rnd)
                	{
                  	n(i,j)++;
                    goto finEval;
                    }


				//
				// Mortalidad de Adultos
				//
				eventRate +=  p.da*N(i,j)/totalRate;
                if(eventRate>=rnd)
                	{
                  	N(i,j)--;
                    goto finEval;
                    }

           	}
		finEval:;
        rnd = ranf();

        time+= -log(rnd)/totalRate;
		dTime = time;
		

	}
	return 0;
}


double DeterministicInvasion(parameters p,simplmat <long> &nTotal,simplmat <long> &NTotal,
				simplmat <long> &nOcupado,simplmat <long> &NOcupado)
{
	double rnd1,rnd2;
	int i,j,priVez=1;

	simplmat <double> N(p.DimX,p.DimY);
	simplmat <double> n(p.DimX,p.DimY);
	
	N.fill(0);
	n.fill(0);
	
#ifdef GRAPHICS
	if(p.gr=='S')
		{
		IGraph(p.DimX,p.DimY);
		ReadIdrisiPalette("Idris256.smp");
		}
#endif
	//
	// Condiciones iniciales fila de 4 celdas con 10 adultos sobre último X
	//
    N(p.DimX-1,p.DimY/2)=p.densIni;
    N(p.DimX-1,(p.DimY/2)-1)=p.densIni;
    N(p.DimX-1,(p.DimY/2)-2)=p.densIni;
    N(p.DimX-1,(p.DimY/2)+1)=p.densIni;

	//
	// Condiciones iniciales (10 adultos en las cuatro celdas centrales)
	//
    //N(p.DimX/2,p.DimY/2)=10;
    //N(p.DimX/2,(p.DimY/2)-1)=10;
    //N((p.DimY/2)-1,p.DimY/2)=10;
    //N((p.DimY/2)-1,(p.DimY/2)-1)=10;
    
    double deltaTime=0.01; // Intervalo de integracion si es = 1 equivale a un modelo discreto
    double sigmaJuv=0.01,sigmaAd=0.001; // SD de fluctuacion estocástica

   	long totalJuv = 0,totalAd = 0;
    double time=0,tempJuv,tempAd,z1,z2,time50=0;
    int dTime=0,dTimeAnt=0;
    int jPlus1,iPlus1,iMinus1,jMinus1;
    long juvOcupado=0,adOcupado=0;
	int salida=0;

	while( 1 )
	{
        totalJuv=0;
        totalAd=0;
        tempJuv=0;
        tempAd=0;
        juvOcupado=0;
        adOcupado=0;
		double cobertura=0;

		// Hacer salida 
		//
		if( dTime>dTimeAnt || priVez || salida)
		{
			for( i=0; i<p.DimX; i++)
				for( j=0; j<p.DimY; j++)
				{
					totalJuv += n(i,j);
					totalAd += N(i,j);
					
					// Criterio para cobertura diámetro de Juv 2m
					//							diámetro de Ad 3.8m
					// cobertura de Ad+Juv>20% en ind/ha
					//
					cobertura=(N(i,j)*M_PI*1.9*1.9+n(i,j)*M_PI)/10000;
					if( cobertura>0.2 )
						adOcupado++;

					if( n(i,j)>0.5 )
						juvOcupado++;
				}
			if(adOcupado >= p.DimX*p.DimY/2 && time50==0)
				time50 = time;  	// Guarda el tiempo al llegar al 50% de ocupación
			//	salida=1;

			nTotal(dTime)= totalJuv;
			NTotal(dTime)= totalAd;
			nOcupado(dTime) = juvOcupado;
			NOcupado(dTime) = adOcupado;
			dTimeAnt = dTime;

			//
			// Hacer salida en forma de mapa
			//

            if( p.writeMap>0 && p.writeMap == dTime )
            {
            	RWFile file;
				ostringstream fname;
				fname << p.fName << dTime << ".map"<< ends;
				
				if(!file.WriteSeed(fname.str().c_str(),N, "AD"))
                    cerr << "No pudo escribir salida AD" << endl;
				if(!file.WriteSeed(fname.str().c_str(),n, "JU"))
                    cerr << "No pudo escribir salida JU" << endl;
			}
            if( p.pomac && dTime==p.tFinal)
            {    
                simplmat <double> prom(4);
                ostringstream pname;
                bool priPro=false;
                int vueltas=MuestreoPorCaminataAlAzar(N,n,prom);
                
                //cout << "Vueltas\t AdGrad1\t AdGrad2\t JuvGrad1\t JuvGrad2\n";
                //cout << vueltas << "\t" << prom(0) << "\t" << prom(1) << "\t" << prom(2) << "\t" << prom(3) << endl;


                pname << p.fName << "Pomac.out" << ends;
                
	        	fstream pomacfile;
				pomacfile.open(pname.str().c_str(),ios::in);
				if(!pomacfile )
					priPro=true;
				pomacfile.close();
				pomacfile.clear();

        		pomacfile.open(pname.str().c_str(), ios::out | ios::app);
		        if( !pomacfile )
    			{
				    cerr << "Cannot open file: " << pname << endl;
				    return 0;
	    		}
				if( priPro )
				{
					pomacfile << "g\tdj0\tdj1\talfa\tb\tda\tdensIni\tAdGrad1\tAdGrad2\tJuvGrad1\tJuvGrad2\tAdOcupado" << endl;
					priPro = false;
				}

                pomacfile << p.g << "\t"
    			        << p.dj0 << "\t"
            			<< p.dj1 << "\t"
            			<< p.alfa << "\t"
        				<< p.b << "\t"
        				<< p.da << "\t"
        				<< p.densIni << "\t"
                        << prom(0) << "\t"
                        << prom(1) << "\t"
                        << prom(2) << "\t"
                        << prom(3) << "\t"
                        << adOcupado << endl;
			}   

          	priVez=0;

#ifdef GRAPHICS
			if(p.gr=='S')
	      		for( j=0; j<p.DimY; j++)
    	  			for( i=0; i<p.DimX; i++)
         				PPix(i,j,long(N(i,j))%256);
#endif
			//if(salida)
			//	return time;

			if( dTime==p.tFinal )
				return time50;
		}


			
      	for( i=0; i<p.DimX; i++)
      		for( j=0; j<p.DimY; j++)
			{
				//tempJuv = n(i,j);
				//tempAd  = N(i,j);
             	// Pasa de Juvenil a Adulto
				tempJuv += -p.g*n(i,j);
				tempAd  +=  p.g*n(i,j);
				// Mortalidad de Juveniles
				tempJuv +=  -(p.dj0+p.dj1*N(i,j))*n(i,j);
				// Mortalidad de Adultos
				tempAd  +=  -p.da*N(i,j);

				// Regeneracion local
				tempJuv +=  (1-p.alfa)*p.b*N(i,j);


#ifdef NOWRAP
                iPlus1  = i + 1;
                jPlus1  = j + 1;
                iMinus1 = i - 1;
                jMinus1 = j - 1;
#else
                iPlus1  = (p.DimX + i + 1) % p.DimX;
                jPlus1  = (p.DimY + j + 1) % p.DimY;
                iMinus1 = (p.DimX + i - 1) % p.DimX;
                jMinus1 = (p.DimY + j - 1) % p.DimY;
#endif


#ifndef DELAYREGE
				tempJuv +=  (p.alfa)/8*p.b*(N(iPlus1,jPlus1)+N(iPlus1,j)+N(iPlus1,jMinus1)+
									N(i,jPlus1)+N(i,jMinus1)+
									N(iMinus1,jPlus1)+N(iMinus1,j)+N(iMinus1,jMinus1)
									);
#else
				tempJuv +=  (p.alfa)/8*p.b*(N(iPlus1,jPlus1)+N(iPlus1,j)+N(iPlus1,jMinus1)+
									N(i,jPlus1)+N(i,jMinus1)+
									N(iMinus1,jPlus1)+N(iMinus1,j)+N(iMinus1,jMinus1)
									);
#endif

				// Acá entra la estocacidad en el proceso.
				// 
				//
                
      			rnd1 = ranf();
      			rnd2 = ranf();
      			z1=sqrt(-2*log(rnd1))*cos(2*M_PI*rnd2);
      			z2=sqrt(-2*log(rnd1))*sin(2*M_PI*rnd2);
				
				n(i,j) = (n(i,j)+tempJuv*deltaTime)*exp(z1*sigmaJuv-sigmaJuv*sigmaJuv/2);
				N(i,j) = (N(i,j)+tempAd*deltaTime)*exp(z2*sigmaAd-sigmaAd*sigmaAd/2);
				
				//n(i,j) = pow(sqrt(n(i,j)+tempJuv*deltaTime)+z1*sigmaJuv,2);
				//N(i,j) = pow(sqrt(N(i,j)+tempAd*deltaTime)+z2*sigmaAd,2);
				
				if(n(i,j)<0.01)
           			n(i,j)=0;

				if(N(i,j)<0.01)
           			N(i,j)=0;

				tempJuv=0;
				tempAd=0;

			}

		
		time+= deltaTime;
		dTime = time;
		
	}
	return 0;
}

// Observador virtual que parte desde el foco de la invasion (dimX, dimY/2) 
// muestrea 10 parcelas
// y asigna a categorias de invasion segun un criterio de densidad de adultos
// Si es menor de 700 ind/ha asigna grado 1 (leve)
// Si es mayor de 1100 ind/ha asigna grado 2 (severo)
// retorna el promedio de adultos y juveniles por grado de invasion
//  prome(0) : Adultos grado1
//  prome(1) : Adultos grado2
//  prome(2) : Juveniles grado1
//  prome(3) : Juveniles grado2
//
int MuestreoPorCaminataAlAzar(simplmat <double> &NN,simplmat <double> &nn,
				simplmat <double> &prome)
{
	double rnd1,rnd2;
	long dimY=NN.getCols()-1,dimX=NN.getRows()-1;
	long pasoX=dimX, pasoY=dimY/2, px=0,py=0;
	int cantGrado1=0, cantGrado2=0;
	int dioVuelta=0;
	prome.fill(0.0);
	
	while(1)
	{
        while(1)
        {
            px=ignuin(0,2)-2;
            py=ignuin(0,4)-2;
            if(px!=0 || py!=0) break;
        }
        
		pasoX += px;
		pasoY += py;
		if( pasoX<0 )
        {
            pasoX=dimX;
            dioVuelta++;
        }
        
		if( pasoY>dimY)
		{
			pasoY=dimY/2;
			dioVuelta++;
		}
		else 
			if(pasoY<0) 
			{
				pasoY=dimY/2;
				dioVuelta++;
			}
		
		if( NN(pasoX,pasoY)<700 && NN(pasoX,pasoY)>100 && cantGrado1<10)
		{
			prome(0) +=NN(pasoX,pasoY);
			prome(2) +=nn(pasoX,pasoY);
			cantGrado1++;
            rnd1= NN(pasoX,pasoY);
            rnd1= prome(0);
		}

		if( NN(pasoX,pasoY) >1100 && cantGrado2<10)
		{
			prome(1) +=NN(pasoX,pasoY);
			prome(3) +=nn(pasoX,pasoY);
			cantGrado2++;
            rnd2 = NN(pasoX,pasoY);
            rnd2 = prome(1);
		}

		if( cantGrado2==10 && cantGrado1==10 || dioVuelta>20 )
		{
			// si dio 20 vueltas y no encontró 10 parcelas se fija que la cantidad sea 7 o mas
			// si no es así asigna 0 a los promedios.
			if( dioVuelta >20 )
				if( cantGrado2<7 || cantGrado1<7 )
				{
					prome(0) = 0;
					prome(1) = 0;
					prome(2) = 0;
					prome(3) = 0;
					return dioVuelta;
				}

			prome(0) /= cantGrado1;
			prome(1) /= cantGrado2;
			prome(2) /= cantGrado1;
			prome(3) /= cantGrado2;
			return dioVuelta;
		}

	}
}

	
int ReadPomacParms(ifstream &inLineFile,parameters &p )
{
	string buff;
	static bool privez=true;
	if(privez)
	{
		inLineFile.open("pomac.lin");
		if( !inLineFile )
		{
			cerr << "Cannot open Parms file" << endl;
			return 0;
		}
		privez=false;
		getline(inLineFile,buff);
	}

	getline(inLineFile,buff);
    if( inLineFile.eof() )
        return 0;

	istringstream ins(buff.c_str());
	ins >> p.g;
	ins >> p.dj0;
	ins >> p.dj1;
	ins >> p.alfa;
	ins >> p.b;
	ins >> p.da;
	ins >> p.densIni;
	
    return 1;
}
