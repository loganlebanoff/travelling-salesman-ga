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

double[][] adj;

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public TravelingSalesman () throws java.io.IOException {

		name = "TravelingSalesman";

		Scanner in = new Scanner(new File(Parameters.dataInputFileName));

		adj = new double[Parameters.numGenes][Parameters.numGenes];
		double sumDist = 0;
		for(int i = 0; i < Parameters.numGenes; i ++) {
			for(int j = 0; j < Parameters.numGenes; j ++) {
				double dist = in.nextDouble();
				adj[i][j] = dist;
				if (i > j) {
					sumDist += dist; // Only count each pairwise distance once
				}
			}
		}
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
	
	public double[][] calcSubpathFitnesses(Chromo X) {
		
		// Subpath fitness is the number of edges in the subpath divided by the subpath length, normalized so that
		// all fitness values add up to 1. Higher fitness values are better.

		int[] chromo = X.chromo;
		int L = chromo.length;
		double[][] totDistance = new double[L][L];
		double[][] fitness = new double[L][L];
		double sum = 0;

		// i is first index of subpath, j is last index of subpath
		for (int i = 0; i < L; i ++) {
			for (int j = i + 1; j < L; j ++) {

				if (!(i == 0 && j == L - 1)) { // Don't allow the entire path to be selected
					totDistance[i][j] = totDistance[i][j-1]; // Start with previous subpath (not counting city j)
					totDistance[i][j] += adj[chromo[j-1]][chromo[j]]; // Add city j-1 to j distance

					fitness[i][j] = (j - i) / totDistance[i][j];
					sum += fitness[i][j];
				}
			}
		}
		// Normalize
		for (int i = 0; i < L; i ++) {
			for (int j = i + 1; j < L; j ++) {
				fitness[i][j] /= sum;
			}
		}
		
		return fitness;
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

