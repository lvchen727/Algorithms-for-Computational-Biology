 <<tools`
 (*This is a full implementation of MEME By Chen Lu and Siqi Zhao*)

 siteEMFasta[inputFilename_, motifWidth_, numMotifs_, maxIterations_, accuracy_]:=
   (*Initializaing three inputs*)
 	Module[{ySequences, lengths, subSeqs, siteProbs, siteFreqs,backgroundProbs, backgroundFreqs}, 
 		(* Read in the input file using readInput from tools,
 		   then use the sequences to produce lists of overlapping windows,
 		   each motifWidth long, to pass into siteEM 
 		   PLEASE GROUP THEM so that each sequence creates its own list of ordered overlapping windows
 		   The actual input to siteEM should be a list of lists of lists, eg:
 		   Original seqs:
 		       ABCDEFG
 		       HIJKLMN
 		       OPQRSTU
 		   Input to siteEM w/motifWidth = 4
 		   {
 		   {{A, B, C, D}, {B, C, D, E}, {C, D, E, F}, {D, E, F, G}},
 		   {{H, I, J, K}, {I, J, K, L}, {J, K, L, M}, {K, L, M, N}},
 		   {{O, P, Q, R}, {P, Q, R, S}, {Q, R, S, T}, {R, S, T, U}}
 		   }
 		   *)
 		   
 		(*Read sequences from the fastq file*)
		ySequences=readInput[inputFilename, motifWidth];
 		(*lengths is the list of the length of each sequences*)
 		lengths = Length[ySequences];
 		(*Store the subsequences in by moving through the origincal sequence by one step*)
 		subSeqs = Map[Partition[#,motifWidth,1]&,ySequences[[1;;lengths]]];
 		 
 		(* Call siteEM, then return the value for lambda (siteProb) and 
 		   the value for theta (siteFreqs) *)
       	{siteProbs, backgroundProbs,siteFreqs, backgroundFreqs} = siteEM[subSeqs, numMotifs, maxIterations, accuracy];
 		Return[{siteProbs, siteFreqs}]
 	]

 siteEM[sequences_, numMotifs_, maxIterations_, accuracy_]:=
	Module[{i = 0, iterMotifs= 1,old, new, oldSiteProb, oldBackgroundProb, oldSiteFreqs, oldBackgroundFreqs, 
		    newSiteProb, newBackgroundProb, newSiteFreqs, newBackgroundFreqs, totalMotifs={}, erasingFactor},
		oldSiteFreqs = Table[Normalize[RandomReal[{1/4,1/2},4],Total], {Length[Transpose[sequences[[1]]]]}];
		oldBackgroundFreqs = Normalize[RandomReal[{1/4,1/2},4],Total];
		oldSiteProb = RandomReal[];
		oldBackgroundProb = 1- oldSiteProb;
		old = {oldSiteProb, oldBackgroundProb, oldSiteFreqs, oldBackgroundFreqs};
		erasingFactor = Map[1.0&, sequences,{2}];
		
		(*Start a loop to find the number of motifs asked for*)
		While[iterMotifs <= numMotifs,
			(* Loop here until either maxIterations has been reached or the sum of the absolute values of the changes from one 
		   iteration to the next in all estimated parameters is less than accuracy.
		   On each iteration, call updateProbs, passing in the old values, to set the new values. 
		   Set old values to be the same as the new values.
		   *)   
			While[True,
				i++;
				new = updateProbs[sequences,oldSiteProb, oldBackgroundProb, oldSiteFreqs, oldBackgroundFreqs, erasingFactor];
				{newSiteProb, newBackgroundProb, newSiteFreqs, newBackgroundFreqs} = new;
				If[i>=maxIterations || Total[Abs[new-old], 3] <= accuracy, 
			    	(*After a motif is found, append the motif to the list of motifs to be returned,
			    	update erasing factors, re-randomize the motif probs*)
					Break[];
				];
	        	{newSiteProb, newBackgroundProb, newSiteFreqs, newBackgroundFreqs} = new;
				old = new;			
			];
		totalMotifs = Append[totalMotifs, {newSiteProb, newBackgroundProb, newSiteFreqs, newBackgroundFreqs}];
		erasingFactor = erasingFactor * Map[erasing[#, newSiteProb, newBackgroundProb, newSiteFreqs, newBackgroundFreqs, Length[Transpose[sequences[[1]]]]]&, sequences];
		Print["erasingFactor", erasingFactor];
		oldSiteFreqs = Table[Normalize[RandomReal[{1/4,1/2},4],Total], {Length[Transpose[sequences[[1]]]]}];
		(*
		oldBackgroundFreqs = Normalize[RandomReal[{1/4,1/2},4],Total];
		oldSiteProb = RandomReal[];
		oldBackgroundProb = 1- oldSiteProb;
		*)
        iterMotifs ++;       
      ];
      Return[Transpose[totalMotifs]];

]


updateProbs[sequences_, oldSiteProb_, oldBackgroundProb_, oldSiteFreqs_, oldBackgroundFreqs_, erasingFactors_] := 
 	Module[{sitePosteriors, backgroundPosteriors, expectationSite, expectationBackground, newSiteProb, newBackgroundProb, newSiteFreqs, newBackgroundFreqs, normalizedPosteriors,erasedNormalizedPosteriors},
 		(*Create list of posterior probabilities of a bound site 
 		having been picked for each draw by calling your sitePosteriors,
  		which you should paste in to this file.*)
  		(*Created a list of list, in which the sublist is the posteriors for  one list*)
  		sitePosteriors = Map[sitePosterior[#,oldSiteProb,oldBackgroundProb,oldSiteFreqs,oldBackgroundFreqs]&,sequences,{2}];
  		
  		(*Apply the erasing factors and the normalizations to the posteriors*)
  		(*Since the sum of backagound and site posterior probabilty adds up to 1, we have:*)
  		normalizedPosteriors = Map[normalization[#,Length[Transpose[sequences[[1]]]]]&,sitePosteriors[[1;;Length[sitePosteriors]]]];
  		erasedNormalizedPosteriors = normalizedPosteriors*erasingFactors;
		backgroundPosteriors = 1- erasedNormalizedPosteriors;
  		
  		(*Now use the posteriors to calculate expecatation for backround and sites.*)
		expectationSite = Table[siteFreqForEachSequence[sequences[[i]],erasedNormalizedPosteriors [[i]]], {i, Length[sequences]}] + 1;
		expectationBackground = Table[siteFreqForEachSequence[sequences[[i]],backgroundPosteriors[[i]]], {i, Length[sequences]}] + 1;
		
		(*Finally, use these counts to compute maximum likelihood estimates for the
		parameters and return these estimates in a list.*)
		(*After calculating the expectation for sites and background, the next step is to calculate the updated probabilities
		and frquencies*)
		
		newSiteProb = Total[Flatten[erasedNormalizedPosteriors,1]]/Length[Flatten[sequences,1]];
		newBackgroundProb = Total[Flatten[1-erasedNormalizedPosteriors,1]]/Length[Flatten[sequences,1]];
		
		newSiteFreqs = Map[Normalize[#, Total] &, Total[expectationSite]];
		newBackgroundFreqs = Total[Map[Normalize[#, Total] &, Total[expectationSite]]]/Length[Transpose[sequences[[1]]]];

		{newSiteProb,newBackgroundProb,newSiteFreqs,newBackgroundFreqs}

(* Make sure to include your siteSample and sitePosterior functions here.*)
]


(* Erasing  to allow for finding multiple motifs. For more detailed explanations of how to implement this, refer to section 2.3 of: [Bailey and Elkan, 1993] Timothy L. Bailey and Charles Elkan. Unsupervised learning of multiple motifs in biopolymers using expectation maximization. Technical Report CS93-302, Department of Computer Science, University of California, San Diego, August 1993.*)
(* Renormalization constraint to reduce bias towards single and double nucleotide repeats (eg. AAAAAAAAAA, TATATATATATA). *)

erasing[subSeq_, siteProb_, backgroundProb_, siteFreqs_, backgroundFreqs_, motifWidth_] :=
Module[{sitePosteriors={}, factor ={}, i = 1, p = 1, eraseFactors = {}, temp = {}},
	sitePosteriors = Map[sitePosterior[#,siteProb,backgroundProb,siteFreqs,backgroundFreqs]&,subSeq];
	factor = 1 - sitePosteriors;	 
	While[i <= motifWidth - 1, p = p*factor[[i]];
   		eraseFactors = Flatten[Append[eraseFactors, p], 1];
   		i++;
   		];
  	While[i <= Length[factor] - motifWidth+1, 
  		 eraseFactors = Flatten[Append[eraseFactors, Apply[Times, factor[[i ;; i + motifWidth-1]]]], 1];
   		 i++;
   		 ];
  		 i = Length[factor];
  		 p = 1;
  	While[i >= Length[factor] - motifWidth + 2, p = p*factor[[i]];
  	    temp = Flatten[Append[temp, p], 1];
        i--;
        ];
        eraseFactors = Flatten[Append[eraseFactors, Reverse[temp]], 1];
   Return[eraseFactors]
   ]
  
(*normalize sitePosteriors for each sequence*) 
normalization[sitePosteriors_, motifWidth_] :=
  Module[{updatedSitePosteriors = {}, i = 1},
  	(*Normalize each position with the posteriors that account for the position*)
   	While[i <= Length[sitePosteriors] - 2* motifWidth + 2,             	
    updatedSitePosteriors = Append[updatedSitePosteriors, Normalize[sitePosteriors[[i ;; i + motifWidth]], Total][[1]]];
    i++;
    ];
    (*The last motifWidth postions cannot be normalized by multiple posteriors, so we normalized them with themselves*)
   updatedSitePosteriors = Append[updatedSitePosteriors, Normalize[sitePosteriors[[Length[sitePosteriors] - 2*motifWidth + 3 ;; Length[sitePosteriors]]], Total]];
   Return[Flatten[updatedSitePosteriors,1]];
   ]

(*The function returns a matrix size of 4 * number of site positions for each sequence;
 For each position, the total frequency of ATCG is normalized to 1*)

siteFreqForEachSequence[subSeq_, sitePosteriorForEachSequence_] :=
 Map[baseFreqForEachSitePosition[#1, sitePosteriorForEachSequence] &, Transpose[subSeq]]

(*This function returns the frequency of ATCG for each position in the sequence*)

baseFreqForEachSitePosition[bases_, sitePosteriorForEachSequence_] :=
 Module[{freq = {0, 0, 0, 0}, i}, 
  Do[freq[[bases[[i]]]] = 
    freq[[bases[[i]]]] + sitePosteriorForEachSequence[[i]], {i, 
    Length[bases]}];
  Return[freq]]

sitePosterior[sequence_, sitePrior_, backgroundPrior_, siteProbs_, backgroundProbs_] := 
 	Module[{i},
 		(*This perhaps requires a lot of explnanation... But a line of code looks cool... It is exactly like dicePosterior execpt the fact that siteProb now is a matrix
 		 The difference is as in diceSample that we need to create a matrix for background to match that matrix. The general formula is the same of the form 1/(1+uodateProbSite/updatedProbBackground)
 		 Times[backgroundProbs[[sequence[[i]]]],{i,Length[sequence]}] means the that we multiple the ith position in the sequence with the probability, so as with the siteProb*)
 		 (sitePrior*Product[siteProbs[[i,sequence[[i]]]],{i,Length[sequence]}])/((sitePrior*Product[siteProbs[[i,sequence[[i]]]],{i,Length[sequence]}]) + (backgroundPrior*Product[backgroundProbs[[sequence[[i]]]],{i,Length[sequence]}]))
 	]
 	
siteSample[siteProb_, backgroundProb_, siteFreqs_, backgroundFreqs_, 
  numDraws_] := 
 Module[{bases = {1,2,3,4}, distSeq, distSites, listPicks, i},
		(* Create a list of draws according to their probabilities, each sequence is assigned to a number as choice*)
		distSeq = EmpiricalDistribution[{siteProb, backgroundProb}->{1,2}];
		(*Pick from the distribution numDraws times of sequences*)
		listPicks = Thread[RandomVariate[distSeq, numDraws]];
		(* Create the distributions of sides for each sequence and put them in a matrix*)
		distSites = {Map[EmpiricalDistribution[#-> bases]&, siteFreqs], Table[EmpiricalDistribution[backgroundFreqs -> bases],{Length[siteFreqs]}]};
		(* Create the output that generates a matrix in which each row is a length four sequence, and the row number is the number of draws*)
		Map[Table[RandomVariate[distSites[[#,i]]], {i,Length[siteFreqs]}]&, listPicks]
 ] 
 