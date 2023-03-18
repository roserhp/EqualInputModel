--General 3 leaves
restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
R=K[b_1,c_1,e_1,b_2,c_2,e_2,b_3,c_3,e_3];
i=1;
m=matrix{{0,b_i,c_i,b_i},{b_i,0,b_i,e_i},{c_i,b_i,0,b_i},{b_i,e_i,b_i,0}}
OE={A,C,G,T}
mm=mutableMatrix {toList apply(OE,i->p_(OE_((position(OE,j->j==i))%4))*(matrix m)_{position(OE,j->j==i)})}
for i to 3 do mm_(i,i)=1-sum flatten entries (matrix mm)^{i}
M1=matrix mm
i=2
M2=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})
i=3
M3=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})

--Identity at the leaves
restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
R=K;
OE={A,C,G,T}
M1=id_(R^4)
M2=M1
M3=M1

S=sort elements (set OE)^**3/splice
netList S

pp=mutableMatrix(R,64,1)
for i to 63 do (
	pp_(i,0)=sum apply(OE,k->p_k*M1_(position(OE,l->l==k),position(OE,l->l==(S_i)_0))
	                           *M2_(position(OE,l->l==k),position(OE,l->l==(S_i)_1))
		                   *M3_(position(OE,l->l==k),position(OE,l->l==(S_i)_2)))) 
pp=matrix pp;

------------------------------------------------------------------------------------------
------------------------- ORIGINAL BASIS -------------------------------------------------
------------------------------------------------------------------------------------------
length select(flatten entries pp,i->i==0) --0
length unique(flatten entries pp) --52 --> 64-52=12 equalities

equalEntries=select(unique for i to 63 list S_(positions(flatten entries pp,l->l==pp_(i,0))),k->length k>1);
length equalEntries --12
netList equalEntries
"3leavesOriginalBasisEqualEntries.txt" << netList equalEntries << endl << close
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

--NEW BASIS
H=transpose(matrix{{1,1,1,1},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
H3=(transpose H)**(transpose H)**(transpose H);
pbar=time H3*pp;

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0))
length unique nonZeroEntries --19
zeroEntries=S_(positions(flatten entries pbar,i->i==0))
length zeroEntries --45

"3leavesBasisHzeroEntries.txt" << netList zeroEntries << endl << close

--PARAMETRIZATION ON THE EIGENVALUES
Rgeneral=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T)]
QBAR=sub(pbar^(positions(flatten entries pbar,i->i!=0)),Rgeneral);
toString QBAR
PBAR=toList apply(nonZeroEntries,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*QBAR_(position(nonZeroEntries,j->j==i),0));
netList PBAR

varp=toList apply(nonZeroEntries,i->(symbol p)_i);
S=K[varp]; 
f=map(Rgeneral,S,PBAR);
f(p_(A,A,A))
f(p_(T,T,T))
I=time trim kernel f;
--16 sec aprox
"3leavesVanishingIdeal.txt" << toString I << endl << close

------------------------------------------------------------------------------------------
---------------------------- STUDY OF THE VANISHING IDEAL --------------------------------
------------------------------------------------------------------------------------------
restart
nonZeroEntries={(A,A,A),(A,C,C),(A,G,G),(A,T,T),(C,A,C),(C,C,A),(C,C,C),(C,C,G),(C,G,C),(G,A,G),(G,C,C),(G,G,A),(G,G,G),(G,T,T),(T,A,T),(T,G,T),(T,T,A),(T,T,G),(T,T,T)}
varp=toList apply(nonZeroEntries,i->(symbol p)_i);
K=frac(QQ[p_A,p_C,p_G,p_T]);
S=K[varp]; 
I=value get "3leavesVanishingIdeal.txt";
betti I
netList I_*

--Identity point: lambda's=1 in the parametrization
QBAR=matrix {{p_A+p_C+p_G+p_T}, {(p_C+p_T)/(p_C*p_T)}, {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(p_A+p_G)/(p_A*p_G)}, {(p_C+p_T)/(p_C*p_T)}, {(p_C+p_T)/(p_C*p_T)},
      {(-p_C^2+p_T^2)/(p_C^2*p_T^2)}, {(-1)/(p_C*p_T)}, {(-1)/(p_C*p_T)}, {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(-1)/(p_C*p_T)},
      {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(-p_A^2+p_C^2-2*p_A*p_G-p_G^2+2*p_C*p_T+p_T^2)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C
      *p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2)}, {1/(p_A*p_G)}, {(p_A+p_G)/(p_A*p_G)}, {1/(p_A*p_G)}, {(p_A+p_G)/(p_A*p_G)}, {1/(p_A*p_G)}, {(-p_A^2+p_G^2)/(p_A^2*p_G^2)}};
sub(I,transpose QBAR)

-- deg 2
I2=ideal toList apply(0..8,i->I_i);
netList I2_*
dim I2, codim I2, degree I2
--(13,6,27)
--deg 5
I5=ideal toList apply(38..40,i->I_i);
netList I5_*
dim I5, codim I5, degree I5
--(17,2,4)
I3=ideal toList apply(9..37,i->I_i);
netList I3_*
dim I3, codim I3, degree I3
--(11,8,12)

codim(I2+I3)--9
codim(I3+I5)--8
codim(I2+I5)--7

CI=ideal{I_0,I_1,I_2,I_3,I_4,I_6,I_13,I_28,I_37};
codim CI --9

jac=jacobian I;
QBAR=matrix {{p_A+p_C+p_G+p_T}, {(p_C+p_T)/(p_C*p_T)}, {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(p_A+p_G)/(p_A*p_G)}, {(p_C+p_T)/(p_C*p_T)}, {(p_C+p_T)/(p_C*p_T)},
      {(-p_C^2+p_T^2)/(p_C^2*p_T^2)}, {(-1)/(p_C*p_T)}, {(-1)/(p_C*p_T)}, {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(-1)/(p_C*p_T)},
      {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(-p_A^2+p_C^2-2*p_A*p_G-p_G^2+2*p_C*p_T+p_T^2)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C
      *p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2)}, {1/(p_A*p_G)}, {(p_A+p_G)/(p_A*p_G)}, {1/(p_A*p_G)}, {(p_A+p_G)/(p_A*p_G)}, {1/(p_A*p_G)}, {(-p_A^2+p_G^2)/(p_A^2*p_G^2)}};
jacId=sub(jac,transpose QBAR);
time rank jacId --9
rank jacId_{0,1,2,3,4,5,6,7,8} --6

--Compute all possible regular sequences made of 6 gens of degree 2 and 3 of degree 3
ind3=time subsets(toList(9..37),3);
ind2=time subsets(toList(0..8),6);
codims3=time apply(toList ind3,i->codim ideal apply(i,j->I_j));
codims2=time apply(toList ind2,i->codim ideal apply(i,j->I_j));
c3=ind3_(positions(codims3,i->i==3));
c2=ind2_(positions(codims2,i->i==6));
c=flatten toList apply(c2,i->toList apply(c3,j->join(i,j)));
codims=time apply(c,i->codim ideal apply(i,j->I_j));
--Process M2 exited abnormally with code 137
--Alternative computation (same but separate for each of the 20 possible regular sequences in deg 2)
c=toList apply(c2,i->toList apply(c3,j->join(i,j)));
cod0=time apply(c_0,i->codim ideal apply(i,j->I_j));
-- used 928.037 seconds
regseq0=(c_0)_(positions(cod0,i->i==9));
netList regseq0
--{0,1,2,3,4,6,13,28,37}
cod1=time apply(c_1,i->codim ideal apply(i,j->I_j));
-- used 989.753 seconds
regseq1=(c_1)_(positions(cod1,i->i==9));
netList regseq1
--++
cod2=time apply(c_2,i->codim ideal apply(i,j->I_j));
-- used 1034.05 seconds
regseq2=(c_2)_(positions(cod2,i->i==9));
netList regseq2
cod3=time apply(c_3,i->codim ideal apply(i,j->I_j));
-- used 942.23 seconds
regseq3=(c_3)_(positions(cod3,i->i==9));
netList regseq3
--++
cod4=time apply(c_4,i->codim ideal apply(i,j->I_j));
-- used 1163.72 seconds
regseq4=(c_4)_(positions(cod4,i->i==9));
netList regseq4
--++
cod5=time apply(c_5,i->codim ideal apply(i,j->I_j));
-- used 889.025 seconds
regseq5=(c_5)_(positions(cod5,i->i==9));
netList regseq5
--++
cod6=time apply(c_6,i->codim ideal apply(i,j->I_j));
-- used 1052.01 seconds
regseq6=(c_6)_(positions(cod6,i->i==9));
netList regseq6
--++
cod7=time apply(c_7,i->codim ideal apply(i,j->I_j));
-- used 1081.96 seconds
regseq7=(c_7)_(positions(cod7,i->i==9));
netList regseq7
--++
cod8=time apply(c_8,i->codim ideal apply(i,j->I_j));
-- used 964.74 seconds
regseq8=(c_8)_(positions(cod8,i->i==9));
netList regseq8
--++
cod9=time apply(c_9,i->codim ideal apply(i,j->I_j));
-- used 900.055 seconds
regseq9=(c_9)_(positions(cod9,i->i==9));
netList regseq9
--++
cod10=time apply(c_10,i->codim ideal apply(i,j->I_j));
 -- used 1197.04 seconds
regseq10=(c_10)_(positions(cod10,i->i==9));
netList regseq10
--++
cod11=time apply(c_11,i->codim ideal apply(i,j->I_j));
-- used 1106.49 seconds
regseq11=(c_11)_(positions(cod11,i->i==9));
netList regseq11
--++
cod12=time apply(c_12,i->codim ideal apply(i,j->I_j));
-- used 1109.4 seconds
regseq12=(c_12)_(positions(cod12,i->i==9));
netList regseq12
--++
cod13=time apply(c_13,i->codim ideal apply(i,j->I_j));
-- used 958.066 seconds
regseq13=(c_13)_(positions(cod13,i->i==9));
netList regseq13
--++
cod14=time apply(c_14,i->codim ideal apply(i,j->I_j));
-- used 1052.5 seconds
regseq14=(c_14)_(positions(cod14,i->i==9));
netList regseq14
--++
cod15=time apply(c_15,i->codim ideal apply(i,j->I_j));
-- used 990.426 seconds
regseq15=(c_15)_(positions(cod15,i->i==9));
netList regseq15
--++
cod16=time apply(c_16,i->codim ideal apply(i,j->I_j));
-- used 931.159 seconds
regseq16=(c_16)_(positions(cod16,i->i==9));
netList regseq16
--++
cod17=time apply(c_17,i->codim ideal apply(i,j->I_j));
-- used 1058.15 seconds
regseq17=(c_17)_(positions(cod17,i->i==9));
netList regseq17
--++
cod18=time apply(c_18,i->codim ideal apply(i,j->I_j));
 -- used 965.575 seconds
regseq18=(c_18)_(positions(cod18,i->i==9));
netList regseq18
--++
cod19=time apply(c_19,i->codim ideal apply(i,j->I_j));
-- used 961.429 seconds
regseq19=(c_19)_(positions(cod19,i->i==9));
netList regseq19
--{1,3,5,6,7,8,9,13,37}
CI2=ideal toList apply({1,3,5,6,7,8,9,13,37},i->I_i);
codim CI2
netList CI2_*
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
--Specific pi
Ipi=sub(I,{p_A=>1/2,p_C=>1/3,p_G=>1/8,p_T=>1/24});
Ipi=sub(I,{p_A=>1/5,p_C=>1/3,p_G=>1/7,p_T=>1-p_A-p_C-p_G});
netList Ipi_*
S2=QQ[varp]
Ipi=sub(Ipi,S2);
loadPackage "Depth"
isPrime Ipi
time isCM(S2/Ipi)
time isCM(Ipi)
time depth(S2/Ipi)
time dim Ipi

codim ideal{Ipi_0,Ipi_1,Ipi_2,Ipi_3,Ipi_4,Ipi_6,Ipi_13,Ipi_28,Ipi_37} --9
CIp=ideal{Ipi_0,Ipi_1,Ipi_2,Ipi_3,Ipi_4,Ipi_6,Ipi_13,Ipi_28,Ipi_37};

dec=time minimalPrimes CIp;
length dec
radical(CIp)==CIp

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

