(* Wolfram Language package *)
(* This project is to simulate a system of any number of interacting repressor*)
nSolveRegVect[synthRtVect_, repressionMatrix_, mDegRtConstVect_, translationRtConstVect_, pDegRtConstVect_, initMRNAConcsVect_, initProteinConcsVect_] :=

 Module[{numGenes = Length[synthRtVect], eqns, funs, solns, sRate, i, j},
  (* Use Table to define a list of numGenes equations for m' and numGene equations for p'. Seq 'eqns' to be the list of equations.*)
  (* construct a list of 2n symbolic equations describing the rates of change of concentrations of the mRNA and protein products of n genes.*)
  eqns = Flatten[ {Table[m[i]'[t] == Product[  1/(1 + Transpose[repressionMatrix][[i]][[j]] * p[j][t]), {j, 1, numGenes} ] * synthRtVect [[i]]- mDegRtConstVect[[i]] * m[i][t], {i, 1, numGenes}],
  			       Table[p[i]'[t] == translationRtConstVect[[i]] * m[i][t] - pDegRtConstVect[[i]] * p[i][t], {i, 1, numGenes}],
  			       Table[m[i][0] == initMRNAConcsVect[[i]], {i , 1, numGenes}],
  			       Table[p[i][0] == initProteinConcsVect[[i]], {i , 1, numGenes}]}
  			     ];

  (* Also use Table to list the functions you want NDSolve to solve for. See docs on NDSolve. set 'funs' to be list of functions.*)
  (* construct a list of 2n functions to solve *)
  funs = Flatten[ Transpose[{Table[m[i], {i, 1, numGenes}],
  	               Table[p[i], {i, 1, numGenes}]}]
  	             ];

  solns = NDSolve[eqns, funs, {t, 0, 40}];
  Print[Plot[Evaluate[Flatten[Table[{m[i][t], p[i][t]}, {i, numGenes}], 1]
  	                          /. solns],
  	        {t, 0, 10},
  	        PlotRange -> All]];
  solns]
