//
// Stochastic simulacion of an exotic tree invasion
//
// Se optimizó para hacer la salida tipo "pomac" a de hectareas ocupadas a 13 y 27 años y la salida de densidad a 26 años
// Solamente con el modelo con ruido log normal, el modelo full estocastico se eliminó
// en los parametros se debe poner tFinal 27 años.
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
    double dj2;
    double alfa;
    double alfa1;
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
					datafile << "DimX\tDimY\ttFinal\tg\tdj0\tdj1\talfa\talfa1\tb\tda\tdensIni\tTime\tJuv\tAd\tJuvOcupa\tAdOcupa" << endl;				
					priPro = false;
				}

			for(i=0; i<ap.tFinal+1; i++)
			{
				datafile << ap.DimX << "\t"
					<< ap.DimY << "\t"
					<< ap.tFinal << "\t"
					<< ap.g << "\t"
					<< ap.dj0 << "\t"
					<< ap.dj1 << "\t"
  					<< ap.dj2 << "\t"
					<< ap.alfa << "\t"
					<< ap.alfa1 << "\t"
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
					<< ap.dj2 << "\t"                    
					<< ap.alfa << "\t"
					<< ap.alfa1 << "\t"					
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
    if(strcmp(buff,"dj2")==0)
    	parms >> p.dj2;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"alfa")==0)
    	parms >> p.alfa;
    else
    	return(1);

	parms >> buff;
    if(strcmp(buff,"alfa1")==0)
    	parms >> p.alfa1;
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
	
	simplmat <double> new_N(p.DimX,p.DimY);
	simplmat <double> new_n(p.DimX,p.DimY);

	new_N.fill(0);
	new_n.fill(0);

		
#ifdef GRAPHICS
	if(p.gr=='S')
		{
		IGraph(p.DimX,p.DimY);
		ReadIdrisiPalette("Idris256.smp");
		}
#endif

	//
	// Condiciones iniciales dos focos de invasión último X
	//
    N(p.DimX-1,p.DimY*.4)=p.densIni;
    N(p.DimX-1,(p.DimY*.4)-1)=p.densIni;
    N(p.DimX-1,(p.DimY*.4)+1)=p.densIni;
    N(p.DimX-1,(p.DimY*.6))=p.densIni;

	//
	// Condiciones iniciales fila de 4 celdas con 10 adultos sobre último X
	//
    //N(p.DimX-1,p.DimY/2)=p.densIni;
    //N(p.DimX-1,(p.DimY/2)-1)=p.densIni;
    //N(p.DimX-1,(p.DimY/2)-2)=p.densIni;
    //N(p.DimX-1,(p.DimY/2)+1)=p.densIni;

    
    double deltaTime=0.01; 				// Intervalo de integracion si es = 1 equivale a un modelo discreto
    double sigmaJuv=0.01,sigmaAd=0.001; // SD de fluctuacion estocástica
    
   	long totalJuv = 0,totalAd = 0;
    double time=0,tempJuv,tempAd,z1,z2,time50=0,cobertura;
    int dTime=0,dTimeAnt=0;
    int jPlus1,iPlus1,iMinus1,jMinus1;
    long juvOcupado=0,adOcupado=0,adOcupado13=0,adOcupado23=0;
	int salida=0;
	simplmat <double> prom(4);					// Virtual observer: vector de promedios con resultado de muestreo
	int vueltas;								// Cantidad de vueltas que dio el virtual observer
	
	int neigh=10; 								// Distancia + 1 del entorno para evaluar la funcion potencial inversa
	simplmat <double> migrants(neigh,neigh); 	// Matriz de dispersion
	
	for(i=0; i<neigh; i++)
		for(j=0; j<neigh; j++)
    		{
			migrants(i,j)=1/pow(sqrt(double(i*i+j*j)),p.alfa1);
			}

	double totalMigrants=0;
	for(i=-neigh+1; i<neigh; i++)
		for(j=-neigh+1; j<neigh; j++)
    		{
            if(i!=0 || j!=0)
        		totalMigrants +=migrants(abs(i),abs(j));
        	}
    double factorMigrants = p.alfa / totalMigrants;
	for(i=0; i<neigh; i++)
		for(j=0; j<neigh; j++)
    		{
			migrants(i,j)*=factorMigrants;
        	}
	// Chequea que la suma de migrants debe ser igual a p.alfa ELIMINAR DESPUES DE DEBUGGING
	/*totalMigrants=0;
	for(i=-neigh+1; i<neigh; i++)
    {
		for(j=-neigh+1; j<neigh; j++)
    		{
            if(i!=0 || j!=0)
	        	totalMigrants +=migrants(abs(i),abs(j));
            cout << migrants(abs(i),abs(j)) << "\t";
        	}
        cout << endl;
    }
	*/
    tempJuv=0;
    tempAd=0;
	
	while( 1 )
	{
		// Hacer salida 
		//
		if( dTime>dTimeAnt || priVez || salida)
		{
	        totalJuv=0;
    	    totalAd=0;
        	juvOcupado=0;
    	    adOcupado=0;
			cobertura=0;
    	    
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

//					if( n(i,j)>0.5 )
//						juvOcupado++;
				}
//			if(adOcupado >= p.DimX*p.DimY/2 && time50==0)
//				time50 = time;  	// Guarda el tiempo al llegar al 50% de ocupación
			//	salida=1;

			nTotal(dTime)= totalJuv;
			NTotal(dTime)= totalAd;
//			nOcupado(dTime) = juvOcupado;
			NOcupado(dTime) = adOcupado;
			dTimeAnt = dTime;

			//
			// Hacer salida en forma de mapa
			//
            if( p.writeMap>0 && (dTime %p.writeMap) == 0 )
            {
            	RWFile file;
				ostringstream fname;
				fname << p.fName << dTime << ".map"<< ends;
				
				if(!file.WriteSeed(fname.str().c_str(),N, "AD"))
                    cerr << "No pudo escribir salida AD" << endl;
				if(!file.WriteSeed(fname.str().c_str(),n, "JU"))
                    cerr << "No pudo escribir salida JU" << endl;
			}

			// Guarda valor de ocupacion a 13 años (1984)
			//
			if ( dTime==13 )
				adOcupado13=adOcupado;
			
			// Guarda valor de ocupacion a 23 años (1994)
			//
			if ( dTime==23 )
				adOcupado23=adOcupado;
				
			// Guarda valor de muestreo de densidad a 26 años
			//
			if ( dTime==26 )
                vueltas=MuestreoPorCaminataAlAzar(N,n,prom);
                
			cerr << dTime<< "\t" << totalJuv << "\t" << totalAd	<< "\t" << adOcupado << endl;
			
            if( p.pomac && dTime==p.tFinal)
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
					pomacfile << "g\tdj0\tdj1\tdj2\talfa\talfa1\tb\tda\tdensIni\tAdGrad1\tAdGrad2\tJuvGrad1\tJuvGrad2\tAdOcupado13\tAdOcupado23" << endl;
					priPro = false;
				}

                pomacfile << p.g << "\t"
    			        << p.dj0 << "\t"
            			<< p.dj1 << "\t"
            			<< p.dj2 << "\t"
            			<< p.alfa << "\t"
            			<< p.alfa1 << "\t"            			
        				<< p.b << "\t"
        				<< p.da << "\t"
        				<< p.densIni << "\t"
                        << prom(0) << "\t"
                        << prom(1) << "\t"
                        << prom(2) << "\t"
                        << prom(3) << "\t"
                        << adOcupado13 << "\t"
                        << adOcupado23 << endl;
			}   

          	priVez=0;

#ifdef GRAPHICS
			if(p.gr=='S')
	      		for( j=0; j<p.DimY; j++)
    	  			for( i=0; i<p.DimX; i++)
                    {
       					cobertura=(N(i,j)*M_PI*1.9*1.9+n(i,j)*M_PI)/10000;
					    if( cobertura>0.2 )
                            PPix(i,j,long(cobertura*30)%256);
                    }
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
				// Mortalidad de Juveniles con termino Allee
				tempJuv +=  -(p.dj0+p.dj1*(N(i,j)-p.dj2)*(N(i,j)-p.dj2))*n(i,j);
				// Mortalidad de Adultos
				tempAd  +=  -p.da*N(i,j);

				// Regeneracion local
				tempJuv +=  (1-p.alfa)*p.b*N(i,j);

				int k,kini,kfin,l,lini,lfin;
				// Dispersion 
				//
				kini = i-neigh+1;
                if(kini<0) kini=0;
                kfin = i+neigh;
                if(kfin>p.DimX) kfin=p.DimX;

                lini = j-neigh+1;
                if(lini<0) lini=0;
                lfin = j+neigh;
                if(lfin>p.DimY) lfin=p.DimY;

                for( k=kini; k<kfin; k++)
                	for( l=lini; l<lfin; l++)
                		{
                        if(k!=i || l!=j)
           					tempJuv += p.b*(N(k,l)*migrants(abs(i-k),abs(j-l)));
                    	}


				// Acá entra la estocacidad en el proceso.
				// 
				//
                
      			rnd1 = ranf();
      			rnd2 = ranf();
      			z1=sqrt(-2*log(rnd1))*cos(2*M_PI*rnd2);
      			z2=sqrt(-2*log(rnd1))*sin(2*M_PI*rnd2);
				//n(i,j) = (n(i,j)+tempJuv*deltaTime);
				//N(i,j) = (N(i,j)+tempAd*deltaTime);
				new_n(i,j) = (n(i,j)+tempJuv*deltaTime)*exp(z1*sigmaJuv-sigmaJuv*sigmaJuv/2);
				new_N(i,j) = (N(i,j)+tempAd*deltaTime)*exp(z2*sigmaAd-sigmaAd*sigmaAd/2);
				
				
				if(new_n(i,j)<0.01)
           			new_n(i,j)=0;

				if(new_N(i,j)<0.01)
           			new_N(i,j)=0;

				tempJuv=0;
				tempAd=0;

			}
			
		for( i=0; i<p.DimX; i++)
      		for( j=0; j<p.DimY; j++)
			{
				n(i,j)=new_n(i,j);
				N(i,j)=new_N(i,j);
			}
		
		time+= deltaTime;
		dTime = time;
		
	}
	return 0;
}

// Observador virtual que parte desde el foco de la invasion (dimX, dimY/2) 
// muestrea 10 parcelas sin reposición.
// *** asigna a categorias de invasion segun un criterio de densidad de adultos
// Si es menor de 700 ind/ha asigna grado 1 (leve)
// Si es mayor de 1100 ind/ha asigna grado 2 (severo)
// *** asigna a categorias de invasion segun un criterio de cobertura 
// *** Si es <=30% asigna grado 1 (leve) sino grado 2
//
// retorna el promedio de adultos y juveniles por grado de invasion
//  prome(0) : Adultos grado1
//  prome(1) : Adultos grado2
//  prome(2) : Juveniles grado1
//  prome(3) : Juveniles grado2
//
int MuestreoPorCaminataAlAzar(simplmat <double> &NN,simplmat <double> &nn,
				simplmat <double> &prome)
{
	double rnd1,rnd2,cobertura=0;
	long dimY=NN.getCols()-1,dimX=NN.getRows()-1;
	long pasoX=dimX, pasoY=dimY/2, px=0,py=0;
	int cantGrado1=0, cantGrado2=0;
	int dioVuelta=0,coberturaCero=0;
	simplmat <long> yaMuestreadoG1(10,2);
    simplmat <long> yaMuestreadoG2(10,2);
	bool repitio=false;
	
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
//		if (coberturaCero>10)
//		{
//			pasoX=dimX;
//			dioVuelta++;
//            coberturaCero=0;
//		}
//		else
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
			for(int i=0;i<cantGrado1;i++)
			{
				if( yaMuestreadoG1(i,0)==pasoX && yaMuestreadoG1(i,1)==pasoY)
				{
					repitio=true;
					break;
				}
			}
			if(!repitio)
			{	
				prome(0) +=NN(pasoX,pasoY);
				prome(2) +=nn(pasoX,pasoY);
				yaMuestreadoG1(cantGrado1,0)=pasoX;
				yaMuestreadoG1(cantGrado1,1)=pasoY;
				cantGrado1++;
			}
			else
				repitio=false;
		}

		if( NN(pasoX,pasoY) >1100 && cantGrado2<10)
		{
			for(int i=0;i<cantGrado2;i++)
			{
				if( yaMuestreadoG2(i,0)==pasoX && yaMuestreadoG2(i,1)==pasoY)
				{
					repitio=true;
					break;
				}
			}
			if(!repitio)
			{	
				prome(1) +=NN(pasoX,pasoY);
				prome(3) +=nn(pasoX,pasoY);
				yaMuestreadoG2(cantGrado2,0)=pasoX;
				yaMuestreadoG2(cantGrado2,1)=pasoY;
				cantGrado2++;
			}
			else
				repitio=false;
		}
/*
		// Si la cobertura es mayor al 0.2 y menor a 0.8 grado1
		// Si la cobertura es mayor a 1
		//
		cobertura=(NN(pasoX,pasoY)*M_PI*1.9*1.9+nn(pasoX,pasoY)*M_PI)/10000;
		if( cobertura>0.2 )
		{
			coberturaCero=0;
			if( cobertura<=0.8 )
			{
				prome(0) +=NN(pasoX,pasoY);
				prome(2) +=nn(pasoX,pasoY);
				cantGrado1++;
			}
			else
				if( cobertura>1)
				{
					prome(1) +=NN(pasoX,pasoY);
					prome(3) +=nn(pasoX,pasoY);
					cantGrado2++;
				}
		}
		else coberturaCero++;
*/
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

	
int ReadPomacParms(ifstream &inLineFile, parameters &p )
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
	ins >> p.dj2;
	ins >> p.alfa;
	ins >> p.alfa1;
	ins >> p.b;
	ins >> p.da;
	ins >> p.densIni;
	
    return 1;
}

/*
int ReadPomacParms(parameters &p )
{
	ifstream inLineFile;
	string buff;
	static int privez=0;
	int line=0;
	if(privez==0)
	{
		privez++;
	}
	inLineFile.open("pomac.lin");
	if( !inLineFile )
	{
		cerr << "Cannot open Parms file" << endl;
		return 0;
	}
	getline(inLineFile,buff);

    
    while(line<privez)
    {
		getline(inLineFile,buff);
		if( inLineFile.eof() )
			return 0;
		line++;
	}
	istringstream ins(buff.c_str());
	ins >> p.g;
	ins >> p.dj0;
	ins >> p.dj1;
	ins >> p.dj2;
	ins >> p.alfa;
	ins >> p.b;
	ins >> p.da;
	ins >> p.densIni;
	privez++;
    inLineFile.close();

	
    return 1;
}
*/
