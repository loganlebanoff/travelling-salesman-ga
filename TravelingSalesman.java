/******************************************************************************
*  A Teaching GA					  Developed by Hal Stringer & Annie Wu, UCF
*  Version 2, January 18, 2004
*******************************************************************************/

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.text.*;

public class TravelingSalesman extends FitnessFunction{

/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/


/*******************************************************************************
*                            STATIC VARIABLES                                  *
*******************************************************************************/

double count; 	// Number of subpaths in a chromosome AND number of pairwise distances between cities
double avgDist; // Avg distance between cities
double[][] adj;

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public TravelingSalesman () throws java.io.IOException {

		name = "TravelingSalesman";

		Scanner in;
		String home = Paths.get(System.getProperty("user.home")).toString();
		if (home.equals("C:\\Users\\Kevin")) {
			in = new Scanner(new File("C:\\Users\\Kevin\\Dropbox\\CAP 5512\\EvoHW4\\TSP\\src\\TSP.data"));
		}
		else {
			in = new Scanner(new File(Parameters.dataInputFileName));
		}

		adj = new double[Parameters.numGenes][Parameters.numGenes];
		double sumDist = 0;
		for(int i = 0; i < Parameters.numGenes; i ++) {
			for(int j = 0; j < Parameters.numGenes; j ++) {
				double dist = in.nextDouble();
				adj[i][j] = dist;
				sumDist += dist;
			}
		}
		count = Parameters.numGenes*(Parameters.numGenes - 1) / 2;
		avgDist = sumDist / count; // TODO: This assumes dist(i,j) = dist(j,i) !!
	}

/*******************************************************************************
*                                MEMBER METHODS                                *
*******************************************************************************/

//  COMPUTE A CHROMOSOME'S RAW FITNESS *************************************

	public void doRawFitness(Chromo X) {

		int[] cuts = new int[2];
		cuts[0] = 0;
		cuts[1] = X.chromo.length - 1;
		X.rawFitness = calcRawFitness(X.chromo, cuts);
	}

	public double calcRawFitness(int[] chromo, int[] cuts) {

		double F = 0;
		for (int i = cuts[0]; i < cuts[1]; i++) {
			F += adj[chromo[i]][chromo[i+1]];
		}
		return F;
	}

//  PRINT OUT AN INDIVIDUAL GENE TO THE SUMMARY FILE *********************************

	public void doPrintGenes(Chromo X, FileWriter output) throws java.io.IOException{
		/*
		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getGeneAlpha(i),11,output);
		}
		output.write("   RawFitness");
		output.write("\n        ");
		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getIntGeneValue(i),11,output);
		}
		Hwrite.right((int) X.rawFitness,13,output);
		output.write("\n\n");
		return;
		*/
	}

/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

}   // End of NumberMatch.java *************************************************

