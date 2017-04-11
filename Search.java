/******************************************************************************
*  A Teaching GA					  Developed by Hal Stringer & Annie Wu, UCF
*  Version 2, January 18, 2004
*******************************************************************************/

import java.io.*;
import java.nio.file.Paths;
import java.util.*;

import java.text.*;

public class Search {

/*******************************************************************************
*                           INSTANCE VARIABLES                                 *
*******************************************************************************/

/*******************************************************************************
*                           STATIC VARIABLES                                   *
*******************************************************************************/

	public static ArrayList<Double> runBests;
	public static ArrayList<Double> runAverages;

	public static FitnessFunction problem;

	public static Chromo[] member;
	public static Chromo[] child;

	public static Chromo bestOfGenChromo;
	public static int bestOfGenR;
	public static int bestOfGenG;
	public static Chromo bestOfRunChromo;
	public static int bestOfRunR;
	public static int bestOfRunG;
	public static Chromo bestOverAllChromo;
	public static int bestOverAllR;
	public static int bestOverAllG;

	public static double sumRawFitness;
	public static double sumRawFitness2;	// sum of squares of fitness
	public static double sumSclFitness;
	public static double sumProFitness;
	public static double defaultBest;
	public static double defaultWorst;

	public static double averageRawFitness;
	public static double stdevRawFitness;

	public static int G;
	public static int R;
	public static Random r = new Random();
	private static double randnum;

	private static int memberIndex[];
	private static double memberFitness[];
	private static int TmemberIndex;
	private static double TmemberFitness;

	private static double fitnessStats[][];  // 0=Avg, 1=Best

	public static double[] avgBBLength;
	public static double[] countBBLength;
	
	public static int extinctionFrequency = 50;
	public static double extinctionProbability = 0;
	public static int extinctionType = 0;
	public static boolean useBBFitnessForXover = false;
	public static double probUseBBFitnessForDispMutation = 1;

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/


/*******************************************************************************
*                             MEMBER METHODS                                   *
*******************************************************************************/


/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/
	
	public static void addBBLength(int L) {
		avgBBLength[G] += L;
		countBBLength[G]++;
	}

	public static void main(String[] args) throws IOException{

		Calendar dateAndTime = Calendar.getInstance(); 
		Date startTime = dateAndTime.getTime();
		
		PrintWriter out = new PrintWriter(new File("output33.txt"));
		PrintWriter BBout = new PrintWriter(new File("BBoutput33.txt"));
		

	//  Read Parameter File
//		System.out.println("\nParameter File Name is: " + args[0] + "\n");
		Parameters parmValues = new Parameters(args[0]);

	//  Write Parameters To Summary Output File
		String summaryFileName = Parameters.expID + "_summary2.txt";
		FileWriter summaryOutput = new FileWriter(summaryFileName);
		parmValues.outputParameters(summaryOutput);

	//	Set up Fitness Statistics matrix
		fitnessStats = new double[2][Parameters.generations];
		for (int i=0; i<Parameters.generations; i++){
			fitnessStats[0][i] = 0;
			fitnessStats[1][i] = 0;
		}
		
		avgBBLength = new double[Parameters.generations];
		countBBLength = new double[Parameters.generations];

	//	Problem Specific Setup - For new new fitness function problems, create
	//	the appropriate class file (extending FitnessFunction.java) and add
	//	an else_if block below to instantiate the problem.
 
		if (Parameters.problemType.equals("TSP")){
				problem = new TravelingSalesman();
		}
	
		else System.out.println("Invalid Problem Type");
		
		// TODO: Remove test functions later
//		Search.testSubpathFitnesses();
//		Search.testSwapMutation();
//		Search.testDisplacementMutation(true);
//		Search.testXover();
		
//		System.out.println(problem.name);

	//	Initialize RNG, array sizes and other objects
		r.setSeed(Parameters.seed);
		memberIndex = new int[Parameters.popSize];
		memberFitness = new double[Parameters.popSize];
		member = new Chromo[Parameters.popSize];
		child = new Chromo[Parameters.popSize];
		bestOfGenChromo = new Chromo();
		bestOfRunChromo = new Chromo();
		bestOverAllChromo = new Chromo();

		if (Parameters.minORmax.equals("max")){
			defaultBest = 0;
			defaultWorst = 999999999999999999999.0;
		}
		else{
			defaultBest = 999999999999999999999.0;
			defaultWorst = 0;
		}

		bestOverAllChromo.rawFitness = defaultBest;

//		System.out.println("Run/tGen/tBest/tAvg/tStDev/n");

		double overallBest = 0;
		double overallAverageBest = 0;
		double overallBestStdv = 0;
		runBests = new ArrayList<Double>();
		runAverages = new ArrayList<Double>();
		ArrayList<Integer> optimalIndividualGens = new ArrayList<Integer>();

//		System.out.println(Parameters.numRuns + " " + Parameters.generations);
		double[][] genBests = new double[Parameters.numRuns + 1][Parameters.generations];
		double[][] genAvgs = new double[Parameters.numRuns + 1][Parameters.generations];

		//  Start program for multiple runs
		for (R = 1; R <= Parameters.numRuns; R++){

			bestOfRunChromo.rawFitness = defaultBest;
//			System.out.println();

			//	Initialize First Generation
			for (int i=0; i<Parameters.popSize; i++){
				member[i] = new Chromo();
				child[i] = new Chromo();
			}

			//	Begin Each Run
			double sum3 = 0;
//			System.out.println("Run\tGen\tBest\tAvg\tStdDev");
			boolean foundOpt = false;
			for (G=0; G<Parameters.generations; G++){

				sumProFitness = 0;
				sumSclFitness = 0;
				sumRawFitness = 0;
				sumRawFitness2 = 0;
				bestOfGenChromo.rawFitness = defaultBest;
				
				if (G % extinctionFrequency == 0) {
					// Probability of performing extinction this time around
					if(Search.r.nextDouble() < extinctionProbability) {
						performMassExtinction();
					}
				}

				//	Test Fitness of Each Member
				for (int i=0; i<Parameters.popSize; i++){

					member[i].rawFitness = 0;
					member[i].sclFitness = 0;
					member[i].proFitness = 0;

					problem.doRawFitness(member[i]);

					sumRawFitness = sumRawFitness + member[i].rawFitness;
					sumRawFitness2 = sumRawFitness2 +
						member[i].rawFitness * member[i].rawFitness;

					if (Parameters.minORmax.equals("max")){
						if (member[i].rawFitness > bestOfGenChromo.rawFitness){
							Chromo.copyB2A(bestOfGenChromo, member[i]);
							bestOfGenR = R;
							bestOfGenG = G;
						}
						if (member[i].rawFitness > bestOfRunChromo.rawFitness){
							Chromo.copyB2A(bestOfRunChromo, member[i]);
							bestOfRunR = R;
							bestOfRunG = G;
						}
						if (member[i].rawFitness > bestOverAllChromo.rawFitness){
							Chromo.copyB2A(bestOverAllChromo, member[i]);
							bestOverAllR = R;
							bestOverAllG = G;
						}
					}
					else {
						if (member[i].rawFitness < bestOfGenChromo.rawFitness){
							Chromo.copyB2A(bestOfGenChromo, member[i]);
							bestOfGenR = R;
							bestOfGenG = G;
						}
						if (member[i].rawFitness < bestOfRunChromo.rawFitness){
							Chromo.copyB2A(bestOfRunChromo, member[i]);
							bestOfRunR = R;
							bestOfRunG = G;
						}
						if (member[i].rawFitness < bestOverAllChromo.rawFitness){
							Chromo.copyB2A(bestOverAllChromo, member[i]);
							bestOverAllR = R;
							bestOverAllG = G;
						}
					}
				}

				// Accumulate fitness statistics
				fitnessStats[0][G] += sumRawFitness / Parameters.popSize;
				fitnessStats[1][G] += bestOfGenChromo.rawFitness;

				averageRawFitness = sumRawFitness / Parameters.popSize;
				stdevRawFitness = Math.sqrt(
							Math.abs(sumRawFitness2 - 
							sumRawFitness*sumRawFitness/Parameters.popSize)
							/
							(Parameters.popSize-1)
							);

				sum3 += bestOfGenChromo.rawFitness;
				//print stats for each generation
				genBests[R][G] = bestOfGenChromo.rawFitness;
				genAvgs[R][G] = averageRawFitness;
				out.println(R + "\t" + G + "\t" + (int)bestOfGenChromo.rawFitness + "\t " + averageRawFitness  + "\t" + stdevRawFitness);

				if(!foundOpt && (int)bestOfGenChromo.rawFitness == 200)
				{
					optimalIndividualGens.add(G);
					foundOpt = true;
				}

				// Output generation statistics to screen
//				System.out.println(R + "\t" + G +  "\t" + (int)bestOfGenChromo.rawFitness + "\t" + averageRawFitness + "\t" + stdevRawFitness);

				// Output generation statistics to summary file
				summaryOutput.write(" R ");
				Hwrite.right(R, 3, summaryOutput);
				summaryOutput.write(" G ");
				Hwrite.right(G, 3, summaryOutput);
				Hwrite.right((int)bestOfGenChromo.rawFitness, 7, summaryOutput);
				Hwrite.right(averageRawFitness, 11, 3, summaryOutput);
				Hwrite.right(stdevRawFitness, 11, 3, summaryOutput);
				summaryOutput.write("\n");


		// *********************************************************************
		// **************** SCALE FITNESS OF EACH MEMBER AND SUM ***************
		// *********************************************************************

				switch(Parameters.scaleType){

				case 0:     // No change to raw fitness
					for (int i=0; i<Parameters.popSize; i++){
						member[i].sclFitness = member[i].rawFitness + .000001;
						sumSclFitness += member[i].sclFitness;
					}
					break;

				case 1:     // Fitness not scaled.  Only inverted.
					for (int i=0; i<Parameters.popSize; i++){
						member[i].sclFitness = 1/(member[i].rawFitness + .000001);
						sumSclFitness += member[i].sclFitness;
					}
					break;

				case 2:     // Fitness scaled by Rank (Maximizing fitness)

					//  Copy genetic data to temp array
					for (int i=0; i<Parameters.popSize; i++){
						memberIndex[i] = i;
						memberFitness[i] = member[i].rawFitness;
					}
					//  Bubble Sort the array by floating point number
					for (int i=Parameters.popSize-1; i>0; i--){
						for (int j=0; j<i; j++){
							if (memberFitness[j] > memberFitness[j+1]){
								TmemberIndex = memberIndex[j];
								TmemberFitness = memberFitness[j];
								memberIndex[j] = memberIndex[j+1];
								memberFitness[j] = memberFitness[j+1];
								memberIndex[j+1] = TmemberIndex;
								memberFitness[j+1] = TmemberFitness;
							}
						}
					}
					//  Copy ordered array to scale fitness fields
					for (int i=0; i<Parameters.popSize; i++){
						member[memberIndex[i]].sclFitness = i;
						sumSclFitness += member[memberIndex[i]].sclFitness;
					}

					break;

				case 3:     // Fitness scaled by Rank (minimizing fitness)

					//  Copy genetic data to temp array
					for (int i=0; i<Parameters.popSize; i++){
						memberIndex[i] = i;
						memberFitness[i] = member[i].rawFitness;
					}
					//  Bubble Sort the array by floating point number
					for (int i=1; i<Parameters.popSize; i++){
						for (int j=(Parameters.popSize - 1); j>=i; j--){
							if (memberFitness[j-i] < memberFitness[j]){
								TmemberIndex = memberIndex[j-1];
								TmemberFitness = memberFitness[j-1];
								memberIndex[j-1] = memberIndex[j];
								memberFitness[j-1] = memberFitness[j];
								memberIndex[j] = TmemberIndex;
								memberFitness[j] = TmemberFitness;
							}
						}
					}
					//  Copy array order to scale fitness fields
					for (int i=0; i<Parameters.popSize; i++){
						member[memberIndex[i]].sclFitness = i;
						sumSclFitness += member[memberIndex[i]].sclFitness;
					}

					break;

				default:
					System.out.println("ERROR - No scaling method selected");
				}

		// *********************************************************************
		// ****** PROPORTIONALIZE SCALED FITNESS FOR EACH MEMBER AND SUM *******
		// *********************************************************************
				
				for (int i=0; i<Parameters.popSize; i++){
					member[i].proFitness = member[i].sclFitness/sumSclFitness;
					sumProFitness = sumProFitness + member[i].proFitness;
				}

		// *********************************************************************
		// ************ CROSSOVER AND CREATE NEXT GENERATION *******************
		// *********************************************************************

				int parent1 = -1;
				int parent2 = -1;

				//  Assumes always two offspring per mating
				for (int i=0; i<Parameters.popSize; i=i+2){

					//	Select Two Parents
					parent1 = Chromo.selectParent();
					parent2 = parent1;
					while (parent2 == parent1){
						parent2 = Chromo.selectParent();
					}

					//	Crossover Two Parents to Create Two Children
					randnum = r.nextDouble();
					if (randnum < Parameters.xoverRate){
						Chromo.mateParents(parent1, parent2, member[parent1], member[parent2], child[i], child[i+1]);
					}
					else {
						Chromo.mateParents(parent1, member[parent1], child[i]);
						Chromo.mateParents(parent2, member[parent2], child[i+1]);
					}
				} // End Crossover

				//	Mutate Children
				for (int i=0; i<Parameters.popSize; i++){
					child[i].doMutation();
				}

				//	Swap Children with Last Generation
				for (int i=0; i<Parameters.popSize; i++){
					Chromo.copyB2A(member[i], child[i]);
				}
				
			} //  Repeat the above loop for each generation

			//output stats per run
			double avgBest = sum3/Parameters.generations;
			out.println("Run " + R + ":\tBest " + (int)bestOfRunChromo.rawFitness + "\tGenOccured " + bestOfRunG + "\tAvgBest " + avgBest);
			overallBest = Math.max((int)bestOfRunChromo.rawFitness, overallBest);
			overallAverageBest += (int)bestOfRunChromo.rawFitness;
			
			runBests.add(bestOfRunChromo.rawFitness);
			runAverages.add(avgBest);

			Hwrite.left(bestOfRunR, 4, summaryOutput);
			Hwrite.right(bestOfRunG, 4, summaryOutput);

			problem.doPrintGenes(bestOfRunChromo, summaryOutput);

			out.println();
		} //End of a Run

		for (int g = 0; g < Parameters.generations; g++) {
			avgBBLength[g] /= countBBLength[g];
			BBout.printf("%.2f\n", avgBBLength[g]);
		}

		Hwrite.left("B", 8, summaryOutput);

		problem.doPrintGenes(bestOverAllChromo, summaryOutput);


		double sum1 = 0;
		double sum2 = 0;

		//	Output Fitness Statistics matrix
		summaryOutput.write("Gen                 AvgFit              BestFit \n");
		for (int i=0; i<Parameters.generations; i++){
			Hwrite.left(i, 15, summaryOutput);
			Hwrite.left(fitnessStats[0][i]/Parameters.numRuns, 20, 2, summaryOutput);
			Hwrite.left(fitnessStats[1][i]/Parameters.numRuns, 20, 2, summaryOutput);
			summaryOutput.write("\n");

		}

		//output overall stats
		double overallBestSum = overallAverageBest;
		double overallBestSumSquares = 0;
		for(Double e: runBests)
			overallBestSumSquares += e*e;
//		double overallAverageSum = 0;
//		double overallAverageSumSquares = 0;
//		for(Double e: runAverages)
//		{
//			overallAverageSumSquares += e*e;
//			overallAverageSum += e;
//		}
		overallAverageBest /= Parameters.numRuns;

//		double overallAverageAverageStd = Math.sqrt( Math.abs( overallAverageSumSquares - overallAverageSum*overallAverageSum/Parameters.numRuns)/(Parameters.numRuns-1));
		double overallAverageBestStd = Math.sqrt( Math.abs( overallBestSumSquares - overallBestSum*overallBestSum/Parameters.numRuns)/(Parameters.numRuns-1));
//		double overallAverageAverage = overallAverageSum/Parameters.numRuns;
		double overallAverageBest95CI = 1.96*(overallAverageBestStd/Math.sqrt(Parameters.numRuns));
		double ci95lo = overallAverageBest - overallAverageBest95CI;
		double ci95hi = overallAverageBest + overallAverageBest95CI; 
		out.println("Overall Stats:");
		out.println("Generation\tAverageBest\tAverageBestStdv\tAverageAverage\tAverageAverageStdv");
		for(int j = 0; j < Parameters.generations; j++)
		{
			double genBestSum = 0;
			double genBestSumSq = 0;
			double genAvgSum = 0;
			double genAvgSumSq = 0;
			for(int i = 1; i <= Parameters.numRuns; i++)
			{
				genBestSum += genBests[i][j];
				genBestSumSq += genBests[i][j]*genBests[i][j];
				genAvgSum += genAvgs[i][j];
				genAvgSumSq += genAvgs[i][j]*genAvgs[i][j];
			}

			double genAvgBest = genBestSum/Parameters.numRuns;
			double genAvgAvg = genAvgSum/Parameters.numRuns;
			double genBestStdv = Math.sqrt( Math.abs( genBestSumSq - genBestSum*genBestSum/Parameters.numRuns)/(Parameters.numRuns - 1));
			double genAvgStdv = Math.sqrt( Math.abs( genAvgSumSq - genAvgSum*genAvgSum/Parameters.numRuns)/(Parameters.numRuns - 1));
			out.println(j + "\t" + genAvgBest + "\t" + genBestStdv + "\t" + genAvgAvg + "\t" + genAvgStdv);
		}

		out.println("OverallAvgBest " + overallAverageBest + "\tOverallBestStdv " + overallAverageBestStd + "\tOverallAvgBest95CI +/- " + overallAverageBest95CI + "\tOverallAvgBestCI95 Range " + ci95lo + " " + ci95hi);
		int optIndivs = optimalIndividualGens.size();
		out.println("Optimal Individuals: " + optIndivs);
		if(optIndivs!=0)
		{
			double sumOptGens = 0;
			double sumOptGensSq = 0;
			for(int e: optimalIndividualGens)
			{
				sumOptGens += e;
				sumOptGensSq += e*e;
			}
			double avgOptGen = sumOptGens/optIndivs;
			double optGenStdev = Math.sqrt( Math.abs( sumOptGensSq - sumOptGens*sumOptGens/optIndivs)/(optIndivs - 1));
			double optGenCI95 = 1.96*(optGenStdev/Math.sqrt(optIndivs));
			double optGenCI95lo = avgOptGen - optGenCI95;
			double optGenCI95hi = avgOptGen + optGenCI95;
			out.println("AverageOptimalGen " + avgOptGen + "\tOptimalGenStDev " + optGenStdev + "\tOptimalGenCI95 +/- " + optGenCI95 + "\tOptimalGenCI95 Range " + optGenCI95lo + " " + optGenCI95hi);
		}
		out.close();
		BBout.close();

		summaryOutput.write("\n");
		summaryOutput.close();

		System.out.println();
		System.out.println("Start:  " + startTime);
		dateAndTime = Calendar.getInstance(); 
		Date endTime = dateAndTime.getTime();
		System.out.println("End  :  " + endTime);
		
	} // End of Main Class
	
	public static void performMassExtinction() {
		// Random Regeneration method
		if(extinctionType == 0) {
			int[] randomIndices = pickRandomIndices(0.8);
			for (int i = 0; i < randomIndices.length; i++) {
				int index = randomIndices[i];
				member[index] = new Chromo();
				int a = 0;
			}
		}
	}
	
	public static int[] pickRandomIndices(double percent) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < member.length; i++) {
			list.add(i);
		}
		Collections.shuffle(list, r);
		int[] res = new int[(int)(percent * member.length)];
		for (int i = 0; i < res.length; i++) {
			res[i] = list.get(i);
		}
		return res;
	}

	public static void testSwapMutation() {
		System.out.println();
		System.out.println();
		Chromo X = new Chromo();
		System.out.println("Before swap mutation:");
		printChromo(X);
		X.doSwapMutation(Search.r.nextInt(Parameters.numGenes));
		System.out.println("After swap mutation:");
		printChromo(X);
	}
	
	public static void testDisplacementMutation(boolean useBBFitness) {
		System.out.println();
		Chromo X = new Chromo();
		System.out.println("Before displacement mutation:");
		printChromo(X);
		X.doDisplacementMutation(useBBFitness);
		System.out.println("After displacement mutation:");
		printChromo(X);
	}
	
	public static void testXover() {
		
		System.out.println();
//		int[] parent1 = new int[]{0,1,2,3,4};
//		int[] parent2 = new int[]{2,4,0,1,3};
		Chromo parent1 = new Chromo();
		Chromo parent2 = new Chromo();
		
		int L = parent1.chromo.length;
		
		Chromo child1 = new Chromo();
		Chromo child2 = new Chromo();
		Chromo.performOrderedXover(parent1, parent2, child1, child2, false);
		
		System.out.println("Parents");
		printChromo(parent2);
		printChromo(parent1);

		System.out.println("Children");
		printChromo(child1);
		printChromo(child2);
	}
	
	public static void testSubpathFitnesses() {
		
		System.out.println("\n");
		Chromo X = new Chromo();
		printChromo(X);
		double[][] subpathFitnesses = ((TravelingSalesman)problem).calcSubpathFitnesses(X);
		int L = X.chromo.length;
		for (int i=0; i<L; i++) {
			for (int j=0; j<L; j++) {
				System.out.printf("%.2f ", subpathFitnesses[i][j]);
			}
			System.out.println();
		}
		System.out.println();

	}

	public static void printChromo(Chromo X) {
		
		for(int l=0; l<X.chromo.length; l++)
			System.out.print(X.chromo[l] + " ");
		System.out.println();
	}

}   // End of Search.Java ******************************************************

