/******************************************************************************
*  A Teaching GA					  Developed by Hal Stringer & Annie Wu, UCF
*  Version 2, January 18, 2004
*******************************************************************************/

import java.io.*;
import java.util.*;
import java.text.*;

public class Chromo
{
/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/

	public int[] chromo;
	public double rawFitness;
	public double sclFitness;
	public double proFitness;

/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/

	private static double randnum;

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public Chromo(){

		ArrayList<Integer> temp = new ArrayList<Integer>();
		for (int i=0; i<Parameters.numGenes; i++)
			temp.add(i);
		Collections.shuffle(temp, Search.r);
		chromo = new int[Parameters.numGenes];
		for(int i = 0; i < Parameters.numGenes; i ++)
			chromo[i] = temp.get(i);

		this.rawFitness = -1;   //  Fitness not yet evaluated
		this.sclFitness = -1;   //  Fitness not yet scaled
		this.proFitness = -1;   //  Fitness not yet proportionalized
	}


/*******************************************************************************
*                                MEMBER METHODS                                *
*******************************************************************************/

	//  Get Alpha Represenation of a Gene **************************************

	public String getGeneAlpha(int geneID){
	/*
		int start = geneID * Parameters.geneSize;
		int end = (geneID+1) * Parameters.geneSize;
		String geneAlpha = this.chromo.substring(start, end);
		return (geneAlpha);
	*/
	return "";
	}

	//  Get Integer Value of a Gene (Positive or Negative, 2's Compliment) ****

	public int getIntGeneValue(int geneID){
	/*
		String geneAlpha = "";
		int geneValue;
		char geneSign;
		char geneBit;
		geneValue = 0;
		geneAlpha = getGeneAlpha(geneID);
		for (int i=Parameters.geneSize-1; i>=1; i--){
			geneBit = geneAlpha.charAt(i);
			if (geneBit == '1') geneValue = geneValue + (int) Math.pow(2.0, Parameters.geneSize-i-1);
		}
		geneSign = geneAlpha.charAt(0);
		if (geneSign == '1') geneValue = geneValue - (int)Math.pow(2.0, Parameters.geneSize-1);
		return (geneValue);
	*/
		return 0;
	}

	//  Get Integer Value of a Gene (Positive only) ****************************

	public int getPosIntGeneValue(int geneID){
	/*	
		String geneAlpha = "";
		int geneValue;
		char geneBit;
		geneValue = 0;
		geneAlpha = getGeneAlpha(geneID);
		for (int i=Parameters.geneSize-1; i>=0; i--){
			geneBit = geneAlpha.charAt(i);
			if (geneBit == '1') geneValue = geneValue + (int) Math.pow(2.0, Parameters.geneSize-i-1);
		}
		return (geneValue);
	*/
	return 0;
	}

	//  Mutate a Chromosome Based on Mutation Type *****************************

	public void doMutation(){

		switch (Parameters.mutationType){

		case 1:    
			
			for (int j=0; j<(Parameters.geneSize * Parameters.numGenes); j++){
				double randnum = Search.r.nextDouble();
				if (randnum < Parameters.mutationRateSwap) {
					this.doSwapMutation(j);
				}
			}
			double randnum = Search.r.nextDouble();
			if (randnum < Parameters.mutationRateDisp) {
				randnum = Search.r.nextDouble();
				if (randnum < Search.probUseBBFitnessForDispMutation) {
					this.doDisplacementMutation(true);
				}
				else {
					this.doDisplacementMutation(false);
				}
			}
			break;
		default:
			System.out.println("ERROR - No mutation method selected");
		}
	}
	
	public void doSwapMutation(int pos1) {
		
		int pos2 = pos1;
		while (pos2 == pos1) {
			pos2 = Search.r.nextInt(chromo.length);
		}
		
		int temp = chromo[pos1];
		chromo[pos1] = chromo[pos2];
		chromo[pos2] = temp;
	}

	public void doDisplacementMutation(boolean useBBFitness) {
		
		int L = chromo.length;
		int[] cuts;
		if (useBBFitness) {
			cuts = this.pickSubtourCutsUseBBFitness();
			Search.addBBLength(cuts[1] - cuts[0] + 1);
		}
		else {
			cuts = Chromo.pickSubtourCutsRandom(L);
		}
//		System.out.println("cuts = " + cuts[0] + ", " + cuts[1]);

		int subpathSize = cuts[1] - cuts[0] + 1;
		int newStart = cuts[0];
		while (newStart == cuts[0]) {
			newStart = Search.r.nextInt(L - subpathSize + 1);
		}
//		System.out.println("newStart = " + newStart);
		
		int[] newChromo = new int[L];
		
		int oldIdx = 0;
		int newIdx = 0;
		while (newIdx < L && oldIdx < L) {
			if (oldIdx == cuts[0]) {
				oldIdx += subpathSize;
			}
			if (newIdx == newStart) {
				newIdx += subpathSize;
			}
			newChromo[newIdx] = chromo[oldIdx];
			oldIdx++;
			newIdx++;
		}
		
		for (int i = 0; i < subpathSize; i++) {
			newChromo[newStart + i] = chromo[cuts[0] + i];
		}
		
		chromo = newChromo;
	}
	
	public int[] pickSubtourCutsUseBBFitness() {
		
		int L = this.chromo.length;
		double rand = Search.r.nextDouble();
		double rWheel = 0;
		int[] cuts = new int[2];
		
		double[][] bbFits = ((TravelingSalesman)Search.problem).calcSubpathFitnesses(this);
		for (int i = 0; i < L; i ++) {
			for (int j = i + 1; j < L; j ++) {
				rWheel = rWheel + bbFits[i][j];
				if (rand < rWheel) {
					cuts[0] = i;
					cuts[1] = j;
					return cuts;
				}
			}
		}
		return cuts;
	}

/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

	public static int[] pickSubtourCutsRandom(int L) {

		int cut1 = Search.r.nextInt(L);
		int cut2 = cut1;
		// Cuts must not be the same and must not encompass the entire chromo
		while (cut2 == cut1 || Math.abs(cut2 - cut1) == L - 1) {
			cut2 = Search.r.nextInt(L);
		}
		
		int[] cuts = new int[2];
		cuts[0] = Math.min(cut1, cut2);
		cuts[1] = Math.max(cut1, cut2);
		return cuts;
	}

	//  Select a parent for crossover ******************************************

	public static int selectParent(){

		double rWheel = 0;
		int j = 0;
		int k = 0;

		switch (Parameters.selectType){

		case 1:     // Proportional Selection
			randnum = Search.r.nextDouble();
			for (j=0; j<Parameters.popSize; j++){
				rWheel = rWheel + Search.member[j].proFitness;
				if (randnum < rWheel) return(j);
			}
			break;

		case 3:     // Random Selection
			randnum = Search.r.nextDouble();
			j = (int) (randnum * Parameters.popSize);
			return(j);

		case 2:     //  Tournament Selection
			double pVal = 0.9;
			int c1 = (int)(Search.r.nextDouble() * Parameters.popSize);
			int c2 = (int)(Search.r.nextDouble() * Parameters.popSize);
			if(Search.member[c1].rawFitness < Search.member[c2].rawFitness)
			{
				int swap = c1;
				c1 = c2;
				c2 = swap;
			}
			return Search.r.nextDouble() <= pVal ? c1 : c2;				

		case 4:		//	Rank Selection
			randnum = Search.r.nextDouble();
			for (j=0; j<Parameters.popSize; j++){
				rWheel = rWheel + Search.member[j].proFitness;
				if (randnum < rWheel) return(j);
			}
			break;

		default:
			System.out.println("ERROR - No selection method selected");
		}
	return(-1);
	}

	//  Produce a new child from two parents  **********************************

	public static void mateParents(int pnum1, int pnum2, Chromo parent1, Chromo parent2, Chromo child1, Chromo child2){


		switch (Parameters.xoverType){

		case 1:    
			//TODO: Give xover chance of using BB-fitness aware method (change false to true with some probability)
			if (Search.useBBFitnessForXover) {
				boolean useBBFitness;
				double randnum = Search.r.nextDouble();
				useBBFitness = randnum < 0.5;
				Chromo.performOrderedXover(parent1, parent2, child1, child2, useBBFitness);
			}
			else {
				Chromo.performOrderedXover(parent1, parent2, child1, child2, false);
			}
			break;

		case 2:     //  Two Point Crossover
			break;

		case 3:     //  Uniform Crossover
			break;

		default:
			System.out.println("ERROR - Bad crossover method selected");
		}

		//  Set fitness values back to zero
		child1.rawFitness = -1;   //  Fitness not yet evaluated
		child1.sclFitness = -1;   //  Fitness not yet scaled
		child1.proFitness = -1;   //  Fitness not yet proportionalized
		child2.rawFitness = -1;   //  Fitness not yet evaluated
		child2.sclFitness = -1;   //  Fitness not yet scaled
		child2.proFitness = -1;   //  Fitness not yet proportionalized
	}

	//  Produce a new child from a single parent  ******************************

	public static void mateParents(int pnum, Chromo parent, Chromo child){

		//  Create child chromosome from parental material
		child.chromo = Arrays.copyOf(parent.chromo, parent.chromo.length);

		//  Set fitness values back to zero
		child.rawFitness = -1;   //  Fitness not yet evaluated
		child.sclFitness = -1;   //  Fitness not yet scaled
		child.proFitness = -1;   //  Fitness not yet proportionalized
	}
	
	public static void performOrderedXover(Chromo parentChromo0, Chromo parentChromo1, 
			Chromo child0, Chromo child1, boolean useBBFitness) {
		
		int[] parent0 = parentChromo0.chromo;
		int[] parent1 = parentChromo1.chromo;

		int L = parent0.length;
		int[][] children = new int[2][L];
		
		int[] cuts;
		if (useBBFitness) {
			// Give each parent equal chance of choosing a good BB
			double rand = Search.r.nextDouble();
			if (rand < 0.5)
				cuts = parentChromo0.pickSubtourCutsUseBBFitness();
			else
				cuts = parentChromo1.pickSubtourCutsUseBBFitness();
		}
		else {
			cuts = Chromo.pickSubtourCutsRandom(L);
		}
		
//		System.out.println("Cut1 = " + cuts[0]);
//		System.out.println("Cut2 = " + cuts[1]);
		
		boolean[] used0 = new boolean[L];
		boolean[] used1 = new boolean[L];

		// Copy other parent's subpath
		for(int i=cuts[0]; i<=cuts[1]; i++) {
			children[0][i] = parent1[i];
			children[1][i] = parent0[i];
			
			used0[children[0][i]] = true;
			used1[children[1][i]] = true;
		}
		
		int start = (cuts[1] + 1) % L;
		int i = start;
		int j0 = start; 
		int j1 = start;
		while (i != cuts[0]) {

			while (used0[parent0[j0]])
				j0 = (j0 + 1) % L;
			while (used1[parent1[j1]])
				j1 = (j1 + 1) % L;

			children[0][i] = parent0[j0];
			children[1][i] = parent1[j1];
			
			i = (i + 1) % L;
			j0 = (j0 + 1) % L;
			j1 = (j1 + 1) % L;
		}
		
		child0.chromo = children[0];
		child1.chromo = children[1];
	}

	//  Copy one chromosome to another  ***************************************

	public static void copyB2A (Chromo targetA, Chromo sourceB){

		targetA.chromo = Arrays.copyOf(sourceB.chromo, sourceB.chromo.length);

		targetA.rawFitness = sourceB.rawFitness;
		targetA.sclFitness = sourceB.sclFitness;
		targetA.proFitness = sourceB.proFitness;
		return;
	}

}   // End of Chromo.java ******************************************************


