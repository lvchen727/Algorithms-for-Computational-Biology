(* Mathematica Raw Program *)

(*The output should be a list containing a vector of b parameters and a matrix of c parameters. 
The input should be a vector with 1s for genes that are regulators 
and 0s for genes that are not, and a random seed*)
generateRandomParams[regulatorsMask_, seed_]:= 
 Module[{bVector, cMatrix, numGenes}, 
 	SeedRandom[seed];
 	numGenes = Length[regulatorsMask];
 	bVector = RandomReal[{0.01, 1}, numGenes];
 	cMatrix = Map[If[# == 0, Table[0, numGenes], RandomReal[{0, 1}, numGenes]]&, regulatorsMask];
 	cMatrix = ReplacePart[cMatrix, {j_, j_} -> 0.];
 	Return[{bVector, cMatrix}]
 ]

(* genericEqns generates system of equations for n mutually repressing proteins. They are
   generic in that the parameters remain as indexed symbols b[j] and c[j,i], except that
   c[i,i] is replaced by zero, implementing the ban on auto-repression.*)
genericEqns[n_]:= 
Module[{eqnsAll, i , j},
	eqnsAll = Table[ -1. + b[i] * m[i] * Product[ 1 + c[j, i] * m[j], {j, DeleteCases[Range[1, n], i]} ] == 0, 
   	                {i, 1, n}];
  	Return[eqnsAll]        
]


(* mutantEqns takes a genePresenceMask indicating which genes were deleted in the cells 
   in which these mRNA measurements were taken. It generates a set of generic equations 
   and then replaces the mRNA level variables for the deleted genes with zeros, ensuring 
   that any actual measurments made on those deleted genes are ignored and any solver 
   doesn't have to deal with them.
*)
mutantEqns[genePresenceMask_]:= 
Module[{eqnsAll, numGenes, i, j},
	numGenes = Length[genePresenceMask];
	eqnsAll = Table[If[genePresenceMask[[i]] == 1, 
                    -1. + b[i] * m[i] * Product[If[genePresenceMask[[j]] == 1, 1 + c[j, i] * m[j], 1], {j, DeleteCases[Range[1, numGenes], i]}] == 0], 
   	                {i, 1, numGenes}];
	Return[DeleteCases[eqnsAll, Null]]
]

(* mutantExpVars returns a list of variables representing the expression levels of
   genes that remain in the genotype described by genePresenceMask. *)
mutantExpVars[genePresenceMask_]:=
Module[{numGenes, expVars = {}, i},
	numGenes = Length[genePresenceMask];
	(*expVars: expression variables*)
	expVars = Table[If[ genePresenceMask[[i]] == 1, Append[expVars, m[i]]], {i, 1, numGenes} ];
	Return[Flatten[DeleteCases[expVars, Null]]]
]

(* expressionEqns replaces the symbolic parameters of the generic equations with actual values
   from a parameter matrix, returning a final set of equations with fixed parameters but variables
   for mRNA levels. Solving these for the mRNA variables yields a set of expression levels that
   are consistent with the parameters.*)
expressionEqns[eqns_, {bVector_, cMatrix_}] :=
 Module[{numGenes, i, k, j, newEqns},
  numGenes = Length[bVector];
  newEqns = eqns;
  (*replace with bVector, store equations in newEqns*)
  newEqns = newEqns /. Thread[Table[b[i], {i, numGenes}] -> bVector];
  (*replace with cMatrix*)
  For[k = 1, k <= numGenes, k++, 
   newEqns = 
    newEqns /. Thread[Table[c[k, j], {j, numGenes}] -> cMatrix[[k]]] 
   ];
  newEqns
]

(* expressionProfile finds a set of steady state mRNA expression levels that are consistent with
   a set of expression equations and a genotype as indicated by genePresenceMask. The result of 
   FindRoot is a list of replacement rules for expression variables. These are then used to construct
   a vector of length n containing the expression levels for genes that are present and zeros for
   genes that are deleted.*)
expressionProfile[expressionEqns_, genePresenceMask_]:= 
Module[{expVars = {}, root, i},
	(*FindRoot[f,{x,x0}],  set variables for the FindRoot function as expVars*)
	expVars = Table[If[genePresenceMask[[i]] == 1, Append[expVars, {m[i], 0}]], {i, 1, Length[genePresenceMask]}];
	expVars  = DeleteCases[expVars, Null] //. {x_List} :> x ;
	(*calculate expression level*)
	root = FindRoot[expressionEqns, expVars];
	(*replace with expression level if gene is not delete, otherwise set as 0*)
	Table[m[i], {i, 1, Length[genePresenceMask]}] /. root /.m[n_] -> 0
]

(* expressionMatrix takes parameters and a list of genePresenceMasks, one for each strain that we
   want an expression profile for. For each mask, it makes a set of equations and calls
   expressionProfile with those equations and the mask. This returns a list of rules that is
   applied to the expression level variables to produce actual values, which are returned.
   The return value is a matrix in which each row represents the expression levels for one
   strain (as described by one genePresenceMask). *)
expressionMatrix[params_, genePresenceMasks_]:=
	(*This ensures that b, c, m are undefined here and in all functions called from here.*)
Module[{i, numMasks},
	numMasks = Length[genePresenceMasks];
	Table[expressionProfile[expressionEqns[mutantEqns[genePresenceMasks[[i]]], params], genePresenceMasks[[i]]], {i,numMasks}]
	
]
	
	(* inferParameters takes a list of expression profiles and a 
   matched list of genePresenceMasks, one for each strain that we have an expression 
   profile for. For each profile-mask pair, it makes a set of equations with the c parameters
   as unknown variables. It is important the the right-hand-side of each equation be zero.
   A vector of the left-hand-sides of these equations is then dotted with itself to make a
   sum of squared errors, where the error is the deviation from zero. NMminize is then used to 
   find values of the c[i,j] parameters that are consistent with the equations, or as close to
   consistent as possible. We must use NMiminize instead of NSolve because the number of 
   equations may be larger than the number of variables and NSolve won't accept over-determined 
   system of equations *)

(*Part 4: Parameter inference.*)

inferParameters[regulatorsMask_, geneExpressionProfiles_, genePresenceMasks_]:=
	(*Make sure b, c, m are undefined here and in all functions called from here.*)
	Block[{b, c, m},

	 ]

(* parameterEqns replaces the symbolic mRNA expression levels of the generic equations with actual values
   from a gene expression matrix, returning a final set of equations with fixed mRNA levesl but variables
   for parameters. Solving these for the parameter variables yields a set of parameters that
   are consistent with the expression levels.*)

parameterEqns[eqns_, geneExpressionProfile_]:=

parameterSumSquaredError[eqns_]:=
		
(*SECTION Some auxiliary functions to support testing parameter inference.*)

(* Generate a random set of parameters, generate expression profiles for those 
   parameters in wild-type cells and cells with all genes deleted one at a time.
   Then use these expression levels to infer the parameters. Compare the inferred 
   parameters to the ones used to generate the expression levels. If the 
   difference is less than or equal to 10^-4 replace it by zero. Return the resulting matrix of differences. 
   
   It should contain nothing but zeros.*)		 

testParameterInference[n_, seed_]:=
	With[{regulatorsMask=ConstantArray[1, {n}],
		  genePresenceMasks=allSinglesAndWTPresenceMatrix[n]},
		With[{params=generateRandomParams[regulatorsMask, seed]},
			With[{expressionMatrix=expressionMatrix[params, genePresenceMasks]},
				With[{diff=params-inferParameters[regulatorsMask,
					 	                  		  expressionMatrix,
					 	     			          genePresenceMasks]},
				Round[diff, 0.0001]
				]]]]

allSinglesAndWTPresenceMatrix[n_]:=
	Join[ConstantArray[1,{1,n}],
		 ConstantArray[1, {n,n}] - IdentityMatrix[n]]
		 
