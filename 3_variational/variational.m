(* ::Package:: *)

(* ::Title:: *)
(*Minimizing Functionals (RBF approach)*)
(*J (y) = (y(x) - Subscript[f, 0](x))^2+a Subscript[f, 1](|y'(x)|) + b Subscript[f, 2](|y''(x)|)*)


(* ::Section:: *)
(*Defining functions Subscript[f, 1], Subscript[f, 2] , Subscript[f, 3] and weight constants a & b*)


(* ::Input:: *)
(*Clear[x,x1,x2,x3,x4,L,f0,f1,f2,a,b,y,u];*)
(*f0[x_]=Sqrt[1+x^2] ;*)
(*f1[x_]=x;*)
(*f2[x_]=Exp[x];*)
(*a=1.;*)
(*b=1.;*)


(* ::Section:: *)
(*Computing Lagrangian and left side of Euler-Lagrange equation*)


(* ::Input:: *)
(**)
(*Lag[x1_,x2_,x3_,x4_]=(x2-f0[x1])^2+a f1[Sqrt[x3^2]]+b f2[Sqrt[x4^2]];*)
(*L[x_]=Lag[x,y[x],y'[x],y''[x]];*)
(*Ly[x_]=\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(x2\)]\(Lag[x1, x2, x3, x4]\)\)/.x1->x/.x2->y[x]/.x3->y'[x]/.x4->y''[x];*)
(*Lyx[x_]=\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(x3\)]\(Lag[x1, x2, x3, x4]\)\)/.x1->x/.x2->y[x]/.x3->y'[x]/.x4->y''[x];*)
(*Lyxx[x_]=\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(x4\)]\(Lag[x1, x2, x3, x4]\)\)/.x1->x/.x2->y[x]/.x3->y'[x]/.x4->y''[x];*)
(*EL[x_]=Ly[x]-Lyx'[x]+Lyxx''[x]//FullSimplify (*left side of Euler-Lagrange equation*)*)


(* ::Section:: *)
(*Defining RBF to occupy and u(x) = Subscript[\[CapitalSigma]c, j]\[Phi](|x-Subscript[x, j]|)*)


(* ::Input:: *)
(*n=3;(*Number of RBF's*)*)
(*Col=Range[0,1,1/(n-1)];(*Colocation points are linearly spaced*)*)
(*\[Epsilon]=1.0; (* Shape parameter *)*)
(*\[Phi][r_]=Exp[-(\[Epsilon] r)^2];(*RBF function*)*)
(*u[x_]=\!\( *)
(*\*UnderoverscriptBox[\(\[Sum]\), \(j = 1\), \(n\)]\( *)
(*\*SubscriptBox[\(c\), \(j\)] \[Phi][*)
(*\*SqrtBox[\(\((x - Col[\([\)\(j\)\(]\)])\)^2\)]]\)\)*)


(* ::Section:: *)
(*Replacing y(x)->u(x) in EL equation and evaluation on interior points of domain*)


(* ::Input:: *)
(*out1=EL[x]/.y->u;*)
(*Eq=RandomReal[{0,1},n-2];*)
(*For[i=1,i<=n-2,i++,*)
(*Eq[[i]]=out1/.x->Col[[i+1]] //FullSimplify	*)
(*]*)
(*Join[Eq,{u[0]-f0[0],u[1]-f0[1]}]*)
(**)


(* ::Section:: *)
(*Solving the non-linear system of equations: Approach 1, Newton Method*)


(* ::Input:: *)
(**)
(*F=Join[{u[0]-f0[0],u[1]-f0[1]},Eq];*)
(*sol=FindRoot[F,{{Subscript[c, 1],1},{Subscript[c, 2],2},{Subscript[c, 3],3}},Method->{"Newton","UpdateJacobian"->1},MaxIterations->100]*)
(*nsol={Subscript[c, 1],Subscript[c, 2],Subscript[c, 3]}/.sol;(*Numeric solution*)*)
(*Plot[{f0[x],u[x]/. sol},{x,0,1},PlotTheme->"Detailed",ImageSize->Large]*)


(* ::Section:: *)
(*Solving the non-linear system of equations: Approach 2, Dynamic System*)
(*Applying the trick of made coefficients Subscript[c, j] dependant of time -> Subscript[c, j](t)*)


(* ::Subsection:: *)
(*Time dependence trick*)


(* ::Input:: *)
(*Subscript[u, t][x_]=u[x]/.{Subscript[c, 1]->Subscript[c, 1][t],Subscript[c, 2]->Subscript[c, 2][t],Subscript[c, 3]->Subscript[c, 3][t]}*)
(*Subscript[Eq, t]=Eq/.{Subscript[c, 1]->Subscript[c, 1][t],Subscript[c, 2]->Subscript[c, 2][t],Subscript[c, 3]->Subscript[c, 3][t]};*)


(* ::Subsection:: *)
(*Dirty way to choose \[Alpha] , \[Beta] and \[Gamma] to ensure convergence*)


(* ::Input:: *)
(*(*Choosing \[Alpha], \[Beta] and \[Gamma]*)*)
(*Clear[F,\[Alpha],\[Beta],\[Gamma]];*)
(**)
(*(*Redefining vectorial field*)*)
(*Print["Vectorial Field"]*)
(*F={\[Alpha](u[0]-f0[0]),\[Beta](u[1]-f0[1]),\[Gamma] Eq[[1]]}*)
(**)
(*(*Jacobian*) *)
(*DF=D[F,{{Subscript[c, 1],Subscript[c, 2],Subscript[c, 3]}}];*)
(**)
(*(*Eigenvalues at stationary point*)*)
(*EVs=Eigenvalues[DF/.sol];*)
(**)
(*(*Reduce[Max[Re[EVs]]<0 && -5.\[LessEqual]\[Alpha]\[LessEqual]5. && -5.\[LessEqual]\[Beta]\[LessEqual]5. && -5.\[LessEqual]\[Gamma]\[LessEqual]5.,{\[Alpha],\[Beta],\[Gamma]}]*)*)
(**)
(*(*Dirty trick*)*)
(*l=10.;*)
(*While[True,*)
(*\[Alpha]=RandomReal[{-l,l}];*)
(*\[Beta]=RandomReal[{-l,l}];*)
(*\[Gamma]=RandomReal[{-l,l}];*)
(*If[Max[Re[EVs]]<0,Break[]];*)
(*]*)
(*Print["Values of \[Alpha],\[Beta],\[Gamma] where real part of eigenvalues is negative"]*)
(*{\[Alpha],\[Beta],\[Gamma],Max[Re[EVs]]}*)
(**)


(* ::Subsection:: *)
(*Solving the dynamic system \[PartialD]*)
(*\!\(\*UnderscriptBox[\(c\), \(_\)]\)=F( *)
(*\!\(\*UnderscriptBox[\(c\), \(_\)]\))*)


(* ::Input:: *)
(*(*Solving the dynamic system*)*)
(*Subscript[T, F]=10;*)
(*s=NDSolve[{*)
(*\[Alpha](Subscript[u, t][0]-f0[0])==Subscript[c, 1]'[t],*)
(*\[Beta](Subscript[u, t][1]-f0[1])==Subscript[c, 2]'[t],*)
(*\[Gamma] Subscript[Eq, t][[1]]==Subscript[c, 3]'[t],*)
(*Subscript[c, 1][0]==1,*)
(*Subscript[c, 2][0]==2,*)
(*Subscript[c, 3][0]==3},*)
(*{Subscript[c, 1][t],Subscript[c, 2][t],Subscript[c, 3][t]},{t,0,Subscript[T, F]=10},*)
(*Method->"ImplicitRungeKutta"];*)
(**)
(*(*Plot of solutions*)*)
(*Plot[Evaluate[{Subscript[c, 1][t],Subscript[c, 2][t],Subscript[c, 3][t]}/.s],{t,0,Subscript[T, F]},PlotTheme->"Detailed",ImageSize->Large]*)
(*Print["Solution \!\(\*SubscriptBox[\(c\), \(1\)]\)[t],\!\(\*SubscriptBox[\(c\), \(2\)]\)[t],\!\(\*SubscriptBox[\(c\), \(3\)]\)[t] at t->\[Infinity]"]*)
(*Subscript[coef, \[Infinity]]={Subscript[c, 1][t],Subscript[c, 2][t],Subscript[c, 3][t]}/.s/.t->Subscript[T, F]*)
(**)
(**)
(*(*Plot of final solution*)*)
(*Plot[{f0[x],u[x]/.Subscript[c, 1]->Subscript[coef, \[Infinity]][[1]][[1]]/.Subscript[c, 2]->Subscript[coef, \[Infinity]][[1]][[2]]/.Subscript[c, 3]->Subscript[coef, \[Infinity]][[1]][[3]]},{x,0,1},PlotTheme->"Detailed",ImageSize->Large]*)
(**)
(*(*Plot of vector fields*)*)
(*ws=0.5; (*Window size of VectorPlot3D around F(Subscript[c, 0])=0*)*)
(*VectorPlot3D[{u[0]-f0[0],u[1]-f0[1], Eq[[1]]},{Subscript[c, 1],nsol[[1]]-ws,nsol[[1]]+ws},{Subscript[c, 2],nsol[[2]]-ws,nsol[[2]]+ws},{Subscript[c, 3],nsol[[3]]-ws,nsol[[3]]+ws},PlotLabel->"F[c] without factors",ImageSize->Large]*)
(*VectorPlot3D[{\[Alpha](u[0]-f0[0]),\[Beta](u[1]-f0[1]), \[Gamma] Eq[[1]]},{Subscript[c, 1],nsol[[1]]-ws,nsol[[1]]+ws},{Subscript[c, 2],nsol[[2]]-ws,nsol[[2]]+ws},{Subscript[c, 3],nsol[[3]]-ws,nsol[[3]]+ws},PlotLabel->"F[c] with factors",ImageSize->Large]*)



