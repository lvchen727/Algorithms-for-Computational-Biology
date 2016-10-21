(*In this assignment, I implement the Forward-Backward algorithm. *)
 
(* decode
INPUT
 - observationSeq is a list of observations, e.g., {2,3,4,1,2}
 - states is a list of the state names, e.g., {m, h}
 - alphabet is a list of the HMM alphabet, e.g., {1, 2, 3, 4}
 - emissionMatrix is a matrix of dimensions {Length[states], Length[alphabet]}.  
     emissionMatrix[[i,j]] is the probability of emitting letter j from state i, 
     e.g., {{0.4, 0.1, 0.1, 0.4}, {0.05, 0.4, 0.5, 0.05}}
 - transitionMatrix is a matrix of dimensions {Length[states], Length[states]}.
     transitionMatrix[[i, j]] is the probability of transitioning to state j on one transition starting from state i.
     e.g., {{0.99, 0.01}, {0.01, 0.99}}
 - initialStateProbs is a list of dimensions {Length[states]}
     initialStateProbs[[i]] is the prior probability that the state from which the first observations was
     emitted is state i.  
OUTPUT
- stateSeq is a list of dimensions {Length[observationSeq]}.
  stateSeq[[i]] is the ith state in the the most likely sequence of states, given the observations. 
  e.g., {h,h,m,m,m}.
  *)

(* TODO: Remove the comments surround this decode function template, complete the function, and test it.
*)

decode[observationSeq_, {states_, alphabet_, emissionMatrix_, transitionMatrix_, initialStateProbs_}] := 
	Module[{stateSeq = {}, forwardMatrix, backwardMatrix, finalMatrix, numStates, numObservations, i = 2, 
		j = Length[observationSeq], currentState, currentStateIndex, k = 1, temp, w =1}, 
		
		(*create forward matrix*)
		numStates = Length[states];
		numObservations = Length[observationSeq];
		forwardMatrix = ConstantArray[0,{numObservations, numStates}];
		forwardMatrix[[1]] = MapThread[#1*#2 &, {initialStateProbs, Transpose[emissionMatrix][[observationSeq[[1]]]]}];
		While[ i <=  numObservations,
			forwardMatrix[[i]] = MapThread[#1*#2 &, {Map[Total[forwardMatrix[[i-1]] * #]&, transitionMatrix], Transpose[emissionMatrix][[observationSeq[[i]]]]}];
			i++;
		];	
		
		(*create backward matrix*)
		backwardMatrix = ConstantArray[0,{numObservations, numStates}];
		backwardMatrix[[j]]  = ConstantArray[1, numStates];
		j = j- 1;
		While[ j >=1,
			backwardMatrix[[j]] = Map[Total[backwardMatrix[[j+1]] * # * Transpose[emissionMatrix][[observationSeq[[j]]]]]&, transitionMatrix];
			j--
		];
		
		(*multiply forward and backward matrices, and normalize finalMatrix*)
		finalMatrix =  forwardMatrix * backwardMatrix;
		temp = Total[finalMatrix];
		
		While[w <= Length[observationSeq],
			finalMatrix[[w]] = finalMatrix[[w]]/temp;
			w++;
		];
		
		(*do the traceback directly from above finalMatrix*)
		While[ k <= Length[observationSeq],
			currentStateIndex = Position[Log[finalMatrix[[k]]], Max[Log[finalMatrix[[k]]]]][[1]][[1]];
			currentState = states[[currentStateIndex]];
			stateSeq = Append[stateSeq, currentState];
			k++;
		];
		
		(* Return the state sequence *)
		Return[stateSeq]
	];



(* calculateAccuracy takes a state sequence genereted from mixed2.fa and calculates 
the number of correctly labeled states.  Note: this function only computes the
accuracy for the mixed2.fa observations.

INPUT
 - stateSeq is a list of state sequences, e.g., {h,m,h,m,m}
 
 OUTPUT
 - numCorrectStates = [int] number of correcly labeled states.

*)
	
calculateAccuracy[stateSeq_] := 
	Module[{keyStateSequence, numCorrectStates},
	
	keyStateSequence = Flatten[Characters[ToLowerCase[Import["mixed2key.fa"]]]];
	numCorrectStates = Count[MapThread[Equal, {stateSeq, keyStateSequence}], True]
	];

(* readHMM takes a HMM text file and outputs the state names, the alphabet
the transition matrix, and emission matrix of the HMM

INPUT
 - file is a path to an HMM file (see notebook for the format of HMM files). 

  OUTPUT
 - states = [list] list of the state names, e.g., {m, h}
 - alphabet = [list] list of the HMM alphabet, e.g., {1, 2, 3, 4}
 - emissionMatrix = [matrix of size numStates x numAlphabet] the emission matrix.  
     Element eij = the probability of state i emitting letter j., e.g., {{0.4, 0.1, 0.1, 0.4}, {0.05, 0.4, 0.5, 0.05}}
 - transitionMatrix = [matrix of size numStates x numStates] the transition matrix.
     Element tij = the probability of state i transitioning to state j, e.g., {{0.99, 0.01}, {0.01, 0.99}}

*)

(*Note: this is not exactly how I would hav written readHMM stylewise, but it works so I won't change it for now. -MRB *)

readHMM[file_] := 
	Module[{a, numStates, alphabet, numAlphabet, firstStateIndex, lastStateIndex,
		states, firstStateProbIndex, lastStateProbIndex, initialStateProbs, 
		firstEmissionIndex, lastEmissionIndex, emissionList, emissionMatrix,
		firstTransitionIndex, lastTransitionIndex, transitionList, transitionMatrix}, 
		
	a = Import[file, {"Text", "Words"}];
	
	numStates = ToExpression[a[[1]]]; (* Use ToExpression to convert from character to number *)

	alphabet = Characters[a[[2]]];
	numAlphabet = Length[alphabet];

	firstStateIndex = 3;
	lastStateIndex = firstStateIndex + numStates - 1;
	states = a[[firstStateIndex ;; lastStateIndex]];

	firstStateProbIndex = lastStateIndex + 1;
	lastStateProbIndex = firstStateProbIndex + numStates - 1;
	initialStateProbs = ToExpression[a[[firstStateProbIndex ;; lastStateProbIndex]]];

	firstEmissionIndex = lastStateProbIndex + 1;
	lastEmissionIndex = firstEmissionIndex + numStates*numAlphabet - 1;
	emissionList = ToExpression[a[[firstEmissionIndex ;; lastEmissionIndex]]];
	emissionMatrix = Partition[emissionList, numAlphabet];

	firstTransitionIndex = lastEmissionIndex + 1;
	lastTransitionIndex = firstTransitionIndex + numStates*numStates - 1;
	transitionList = ToExpression[a[[firstTransitionIndex ;; lastTransitionIndex]]];
	transitionMatrix = Partition[transitionList, numStates];
	
	{states, alphabet, emissionMatrix, transitionMatrix, initialStateProbs}

];

	
(* readFasta reads a fasta file and outputs the nucleotide sequence converted to numbers
INPUT
- fastaFile is a string representing the path to fasta file

OUTPUT
- input is a list of bases in the file indicated by fastaFile.  
  bases are translated form ACGT to 1234.
  e.g., {1,3,2,4,2}
*)
readFasta[fastaFile_]:=
	Flatten[Map[Characters, Import[fastaFile]] 
		   /. {"A"->1, "C"->2, "G"->3, "T"->4}
		   ]