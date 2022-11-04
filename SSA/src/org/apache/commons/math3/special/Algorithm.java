package org.apache.commons.math3.special;
import java.util.Scanner;
import java.util.Arrays;

import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.util.ContinuedFraction;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;

public class Algorithm {
	
	public void Data(int fasilitas, int lokasi) {
		try (Scanner input = new Scanner(System.in)) {
			System.out.println("Masukkan banyak fasilitas :");
			fasilitas = input.nextInt();
			System.out.println("Masukkan banyak lokasi :");
			lokasi = input.nextInt();
			
			//input data
			int[][] inputdata = new int[fasilitas][lokasi];
			for (int i = 0; i < fasilitas; i++) {
			    for (int j = 0; i < lokasi; i++ ) {
			    inputdata[i][j] = input.nextInt(); //masukkan data inputan
			    }
			}
			
			//bentuk matriksnya
			for (int i =0; i<fasilitas; i++) {
				for(int j = 0; j<lokasi; j++) {
					System.out.print(inputdata[i][j] + "");
				}
				System.out.println();
			}
		}
	}
	
	//inisialisasi parameter
    private static final int numFS = 5;
    private static final int numDimensions = 30;
    private static final double pdp = 0.1;
    private static final double dg = 0.6; 
    private static final double Gc = 1.9;
    private static final int max_iter = 1;
    private static final double lb = 1.5; //lambda
    private static double tupai[][];
    
    
    //membangkitkan populasi tupai terbang awal
    public void Population(int numFS, int fasilitas, int fsu, int fsl) {
    	fsu=0;
    	fsl=1;
    	for (int i =0; i<numFS; i++) {
    		for (int j =0; j<fasilitas; j++) {
    			tupai = new double[fasilitas][numFS];
				//tupai[i][j] =  (fsl + (Math.random() * ((fsu - fsl)))); 
    			tupai[i][j]= (fsl + (Math.random() * ((fsu - fsl))));
    			
    		}
    	}
    }
    
    //pengurutan tupai
    public void UrutFS(int numFS, int fasilitas, int tupais) {
    	for (int i =0; i<numFS; i++) {
    		for (int j =0; j<fasilitas; j++) {
    			int k = 1;
    			for (int m =0; m<fasilitas; m++) {
    				if (tupai[i][j] > tupai[i][m]) {
    	                k = k+1;
    	            }
    			}
    			int[][] urut;
				urut[i][j] = k;
    		}
    			
    	}
    }
    
    //menentukan matriks penempatan
    public void matriks(int numFS, int fasilitas, int urut) {
    	for (int i =0; i<numFS; i++) {
    		for (int j =0; j<fasilitas; j++) {
    			for (int k=0; k<fasilitas; k++) {
    				int[][] x;
					if (urut[i][k] == j) {
    					x[j][k] = 1;
    				} 
    				else {
    					x[j][k] = 0;
    				}
    					
    			}
    			    			
    		}
    			
    	}
    }
    
  //menghitung nilai fungsi tujuan(nilai fitness)
    public void calcfitness(int numFS, int fasilitas, int frekuensi1,int frekuensi2, int jarak, int x) {
    	for (int i =0; i<numFS; i++) {
    		int[] fitness;
			fitness[i] = 0; //nilaifitness
    		for (int j =0; j<fasilitas; j++) {
    			for (int k=0; k<fasilitas; k++) {
    				for (int l=0; l<fasilitas; l++) {
    					for (int m=0; m<fasilitas; m++) {
    						fitness[i] = fitness[i] + (0.5 * frekuensi1[j][l] * jarak[k][m] * 
    							   x[j][k] * x[l][m]) + 0.5 * (frekuensi2[j][l] * jarak[k][m] * 
    							   x[j][k] * x[l][m]));
    					}
    				}	
    					
    			}
    			    			
    		}
    			
    	}
    }
    
    //menentukan lokasi masing-masing tupai terbang
    public void location(int numFS, int fitness) {
    	int[] tampung2 = null;
		for (int i = 0; i<numFS; i++){
			tampung2[i] = i;
    	}
    	for (int i = 0; i<numFS; i++){
    		for (int j =0; j<numFS-1; j++) {
    			int[] tampung = null;
				if (tampung[j] > tampung[j+1]) {
    				int temp = tampung[j];
    				tampung[j] = tampung[j+1];
    				tampung[j+1] = temp;
    				
    				int temp2 = tampung2[j];
    				tampung2[j] = tampung2[j+1];
    				tampung2[j+1] = temp2;
    			}
    		}
    	}
    }
    
    //modifikasi lokasi tupai terbang berada pada pohon acorn
    public void newloc (int numFS, int fasilitas, int tupai, double dg, double Gc, double pdp) {
    	for (int i = 0; i<numFS; i++){
    		for (int j =0; j<fasilitas; j++) {
    			if (i >=2 && i<=4) {
    				double R1 = Math.random();
    				double[][] tupaibaru;
					if (R1 >= pdp) {
    					tupaibaru [i][j] = tupaibaru[i][j] + dg *Gc * (tupai[1][j] - tupai[i][j]);
    				}
    				else {
    					tupaibaru [i][j] = Math.random();
    				}
    			}
    		}
    	}
    	
    }
    
  //modifikasi lokasi tupai terbang berada pada pohon normal
    public void normaltree (int numFS, int fasilitas, int tupai, double dg, double Gc, double pdp) {
    	for (int i = 0; i<numFS; i++){
    		for (int j =0; j<fasilitas; j++) {
    			double k = Math.random();
    			if (k <= 0.5) {
    				double R2 = Math.random();
    				if (R2 >= pdp) {
    					tupaibaru [i][j] = tupaibaru[i][j] + dg *Gc * (tupai[2][j] - tupai[i][j]);
    				}
    				else {
    					tupaibaru [i][j] = Math.random();
    				}
    			}
    			if (k > 0.5) {
    				double R3 = Math.random();
    				if (R3 >= pdp) {
    					tupaibaru [i][j] = tupaibaru[i][j] + dg *Gc * (tupai[1][j] - tupai[i][j]);
    				}
    				else {
    					tupaibaru [i][j] = Math.random();
    				}
    			}
    		}
    	}
    	
    }
    
    //menghitung dan update nilai konstanta musim dan konstanta musim minimum
    public void Calc_update_konst (int fasilitas, int tupaibaru) {
    	int iterasi;
    	int a = -6;
    	double b = max_iter/2.5;
		if (iterasi <= 3) {
    		double sc_1 = 0;
    		for (int i = 0; i<fasilitas; i++){
    			
    			sc_1 = sc_1 + (Math.pow(tupaibaru[i][j] - tupaibaru[1][j]));
    		}
    		double sc = Math.sqrt(sc_1);
    		//10e^-6/365^iterasi/max_iter/2.5
    		double smin = (Math.pow(Math.exp(10),a))/(Math.pow(365, b));
    	}
    	
    }
    
    //memperbarui lokasi tupai terbang menggunakan levy flight
    public void levy (int numFS, int fasilitas, int tupaibaru, int fsl, int fsu, int lb, double num, int den, double sigma_u, double s) {
    	fsu=0;
    	fsl=1;
    	for (int i = 0; i<numFS; i++){
    		for (int j =0; j<fasilitas; j++) {
    			double q = 1/lb;
    			double r = num/den;
    			double z = 1/Gamma;
    			double w = 1+lb;
    			num = Gamma(1+lb)* Math.sin(Math.PI*lb/2); //pembilang
    			den = Gamma((1+lb)/2)*lb*2^((lb-1)/2); //penyebut
    			sigma_u = (Math.pow(r, q));
    			
    			U[j] = Math.random(); //random(0,1)
    			V[j] = Math.normal(); //normal(0.1)
    			
    			s = U/(Math.pow(Math.abs(V), z)); //U/|V|^1/2
    			L[j] = (lb*Gamma*Math.sin(Math.PI*lb/2)/Math.PI)* 1/(Math.pow(s, w)); //nilai levy
    		}
    		if (sc < smin) {
    			tupaibaru[i][j] = fsl+L[j]*(fsu-fsl);
    		}
    		else {
    			tupaibaru[i][j] = tupaibaru[i][j];
    		}
    	}
    }
    
    //update tupai terbang
    public void updateFS (int numFS, int fasilitas, int tupaibaru, int tupai) {
    	for (int i = 0; i<numFS; i++){
    		for (int j =0; j<fasilitas; j++) {
    			tupai[i][j] = tupaibaru[i][j];
    		}
    	}	
    }
    

	private int Gamma(int i) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	
        
    
    
    
}
