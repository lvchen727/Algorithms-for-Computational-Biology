 (*In these assignments, we have a bag containing two types of dice with different probabilities of rolling each number.  Someone selects a die from the bag at random, rolls it a fixed number of times, reports the outcomes, returns it to the bag, and repeats the process.  Here, I write code to run the EM (Expectation Maximization) algorithm to estimate the parameters of the system -- the probability of drawing each die type and the conditional probability of each face given the die type.*)
  

(*The function diceEM takes a sample of die draws and rolls, runs the EM algorithm, and outputs the estimated parameters.*)
 
 diceEM[sample_, maxIterations_, accuracy_]:=
	Module[{numFaces, binCounts, oldType1Prob, oldType2Prob, oldFaceProbs1, oldFaceProbs2, 
		    newType1Prob, newType2Prob, newFaceProbs1, newFaceProbs2, iterationCount, sumOfChanges},
		(* Initialize the local variables here.*)
		numFaces = Max[sample];
		binCounts = sampleToBinCounts[sample, numFaces];
		oldType1Prob = 0.45;
		oldType2Prob = 0.55;
		oldFaceProbs1 = Normalize[RandomReal[{1/numFaces, 2/numFaces}, numFaces], Total];
		oldFaceProbs2 = Normalize[RandomReal[{1/numFaces, 2/numFaces}, numFaces], Total];
		iterationCount = 1; 
		sumOfChanges = Total[oldFaceProbs1] + Total[oldFaceProbs2] + 1;

		(* Loop here until either maxIterations has been reached or the sum of the absolute values of the changes from one 
		   iteration to the next in all estimated parameters is less than accuracy.
		   On each iteration, call updateProbs, passing in the old values, to set the new values. 
		   Then test whether the termination conditions have been met (Return[] breaks a loop returning its argument).
		   Finally, if termination conditions have not been met, set old values to be the same as the new values.
		   *)
		While[(iterationCount <= maxIterations &&  sumOfChanges > accuracy),
			 {newType1Prob, newType2Prob, newFaceProbs1, newFaceProbs2} = updateProbs[binCounts, oldType1Prob, oldType2Prob, oldFaceProbs1, oldFaceProbs2];
			 
			 sumOfChanges = Total[Abs[newFaceProbs1 - oldFaceProbs1]] + Total[Abs[newFaceProbs2 - oldFaceProbs2]] 
			 				+ Abs[oldType1Prob- newType1Prob] + Abs[oldType2Prob- newType2Prob];
			 				
			 {oldType1Prob, oldType2Prob, oldFaceProbs1, oldFaceProbs2} = {newType1Prob, newType2Prob, newFaceProbs1, newFaceProbs2};
			 iterationCount ++;
			
		]; 
		   
		(*At the end, return the estimated parameters with the less likely die first.*)
		If[newType1Prob <= newType2Prob,
		   {newType1Prob, newType2Prob, newFaceProbs1, newFaceProbs2},
		   {newType2Prob, newType1Prob, newFaceProbs2, newFaceProbs1}]
	]

sampleToBinCounts[sample_, numFaces_] :=
	Module[{eachDraw, eachFace},
			Table[Count[eachDraw, eachFace], {eachDraw, sample}, {eachFace, Range[numFaces]}]
	]
		
updateProbs[binCounts_, oldType1Prob_, oldType2Prob_, oldFaceProbs1_, oldFaceProbs2_] :=
	Module[{posteriors, type1Count, type2Count, faceCounts1, faceCounts2},
		(*Create list of posterior probabilities of a Type1 die having been rolled on each draw by calling your dicePosteriors,
		  which you should paste in to this file. *)
		posteriors = Map[dicePosterior[#1,oldType1Prob, oldType2Prob, oldFaceProbs1, oldFaceProbs2]&,binCounts];
		(* Now use the posteriors to calculate EXPECTED counts for each die and each face in the sample.*)
		faceCounts1 = Total[MapThread[#1 * #2 &, {binCounts, posteriors}]];	
		faceCounts2 = Total[MapThread[#1 * #2 &, {binCounts, 1 - posteriors}]];
		type1Count = Total[Length[binCounts] * posteriors];
		type2Count = Total[Length[binCounts] * (1 - posteriors)];
		
		(* Finally, use these counts to compute maximum likelihood estimates for the parameters and return these estimates
		   in a list.*)
		{type1Count /(type1Count + type2Count), type2Count /(type1Count + type2Count) , 
		Normalize[faceCounts1, Total], Normalize[faceCounts2, Total]}
		
	]
c
(* Make sure to include your diceSample and dicePosterior functions here.*)
dicePosterior[binCounts_, type1Prior_, type2Prior_, faceProbs1_, faceProbs2_] := 
	Module[{(*Pr(Die1, binCounts)*)
	    die1Probability = Apply[Times, MapThread[If[#1==0 && #2 == 0, 1, #1^#2]&, {faceProbs1, binCounts}]]* type1Prior, 
	    (*Pr(Die2, binCounts)*)
	    die2Probability = Apply[Times, MapThread[If[#1==0 && #2 == 0, 1, #1^#2]&, {faceProbs2, binCounts}]]* type2Prior},
	
		(*If Pr(Die1, binCounts) (i.e., die1Probability) is 0, the posterior probability of Die 1 is 0, 
	  	otherwise, it equals 1/(1+die2Probability/die1Probability)*)
		If[die1Probability == 0, 0, 1/ (1 + die2Probability/ die1Probability) ]
	]


diceSample[numType1_, numType2_, type1_, type2_, draws_, rollsPerDraw_] :=
	Module[{},
 		Map[If[# == 0, 
 		   RandomVariate[EmpiricalDistribution[type1 -> Range[Length[type1]]], rollsPerDraw], 
 		   RandomVariate[EmpiricalDistribution[type2 -> Range[Length[type2]]], rollsPerDraw]
 	      ]&,
 	    RandomVariate[BernoulliDistribution[(numType2)/(numType1 + numType2)], draws]]
 	]