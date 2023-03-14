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
-- deg 2
I2=ideal toList apply(0..8,i->I_i);
netList I2_*
dim I2, codim I2, degree I2
--(13,6,27)
betti (trim I2)
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6} --6
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_7} --6
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_6,I2_7} --6
codim ideal{I2_0,I2_1,I2_2,I2_4,I2_6,I2_7} --5
codim ideal{I2_0,I2_1,I2_3,I2_4,I2_6,I2_7} --5
codim ideal{I2_0,I2_2,I2_3,I2_4,I2_6,I2_7} --5
codim ideal{I2_1,I2_2,I2_3,I2_4,I2_6,I2_7} --5
loadPackage "Depth"
isRegularSequence({I2_0,I2_1,I2_2,I2_3,I2_4,I2_6})

--deg 5
I5=ideal toList apply(38..40,i->I_i);
netList I5_*
dim I5, codim I5, degree I5
--(17,2,4)
codim ideal{I5_0,I5_1}--2
codim ideal{I5_0,I5_2}--2
codim ideal{I5_1,I5_2}--2

codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I5_1} --7
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I5_2} --7
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_1,I5_2} --7

codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_7,I5_0,I5_1} --6
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_7,I5_0,I5_2} --7
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_7,I5_1,I5_2} --7

codim ideal{I2_0,I2_1,I2_2,I2_3,I2_6,I2_7,I5_0,I5_1} --6
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_6,I2_7,I5_0,I5_2} --6
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_6,I2_7,I5_1,I5_2} --6

I3=ideal toList apply(9..37,i->I_i);
netList I3_*
dim I3, codim I3, degree I3
--(11,8,12)
codim ideal{I3_0,I3_1,I3_2} --2
codim ideal{I3_18,I3_19,I3_20} --3
for i from 0 to 28 do print(i,codim ideal{I3_18,I3_19,I3_20,I3_i})
codim ideal{I3_18,I3_19,I3_20,I3_23}
for i from 0 to 28 do print(i,codim ideal{I3_18,I3_19,I3_20,I3_23,I3_i})
for i from 0 to 28 do print(i,codim ideal{I3_18,I3_19,I3_20,I3_23,I3_8,I3_i})
codim ideal{I3_18,I3_19,I3_20,I3_23,I3_8,I3_1} --6

codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I5_1,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1}--9
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_1,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1}--8
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1}--8
codim ideal{I2_0,I2_1,I2_2,I5_0,I5_1,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1}--7

aux={I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I5_1,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1}
for i to 13 do print(i,codim ideal drop(aux,{i,i}))


codim ideal {I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I5_1,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1}--9
codim ideal {I2_1,I2_2,I2_3,I2_6,I5_0,I5_1,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1}--8
codim ideal {I2_0,I2_1,I2_2,I2_3,I2_6,I5_0,I5_1,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1}--9

aux2={I2_1,I2_2,I2_3,I2_6,I5_0,I5_1,I3_20,I3_23,I3_8,I3_1}
codim ideal aux2

for i from 0 to 8 do print(i,codim ideal{I3_18,I3_19,I3_20,I3_23,I3_8,I3_1,I2_i}) 
for i from 0 to 2 do print(i,codim ideal{I3_18,I3_19,I3_20,I3_23,I3_8,I3_1,I5_i}) 

for i to 28 do print(i,codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I3_i}) 
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I3_20} --8
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I5_1,I3_20} --8
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I3_28} --8
codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I3_20,I3_28} --8

for i to 28 do print(i,codim ideal{I2_0,I2_1,I2_2,I2_3,I2_6,I5_0,I3_20,I3_i}) 
for i to 28 do print(i,codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I5_0,I3_20,I3_i}) 


codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6} --6
codim ideal{I3_18,I3_19,I3_20,I3_23,I3_8,I3_1} --6

apply({I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I3_18,I3_19,I3_20,I3_23,I3_8,I3_1,I5_0,I5_1},i->codim ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6,I5_0,I3_20,i})

codim(I2+I3)--9
codim(I3+I5)--8
codim(I2+I5)--7


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
CIpi=time regSeqInIdeal(Ipi,9);
netList CIpi_*
isHomogeneous CIpi

betti CIpi
CIpi_6
-- it does find regular sequences but they are not gens, one even not homogeneous

codim ideal{Ipi_0,Ipi_1,Ipi_2,Ipi_3,Ipi_4,Ipi_6,Ipi_13,Ipi_28,Ipi_37} --9
CIp=ideal{Ipi_0,Ipi_1,Ipi_2,Ipi_3,Ipi_4,Ipi_6,Ipi_13,Ipi_28,Ipi_37};

dec=time minimalPrimes CIp;
length dec
radical(CIp)==CIp

jacIpi=jacobian Ipi;
--rank jacIpi
jacIpiId=sub(jacIpi,{p_(A,A,A)=>p_A, p_(A,C,C)=>0, p_(A,G,G)=>0, p_(A,T,T)=>0, p_(C,A,C)=>0, p_(C,C,A)=>0, p_(C,C,C)=>p_C,
       p_(C,C,G)=>0, p_(C,G,C)=>0, p_(G,A,G)=>0, p_(G,C,C)=>0, p_(G,G,A)=>0, p_(G,G,G)=>p_G, p_(G,T,T)=>0,
       p_(T,A,T)=>0, p_(T,G,T)=>0, p_(T,T,A)=>0, p_(T,T,G)=>0, p_(T,T,T)=>p_T})
rank jacIpiId_{0,1,2,3,4,5,6,7,8} --6
rank jacIpiId_{0,2,4,5,7,8} --6
codim ideal{Ipi_0,Ipi_2,Ipi_4,Ipi_5,Ipi_7,Ipi_8}--4
rank jacIpiId --19

jacIpiId=sub(jacIpi,sub(transpose QBAR,{p_A=>1/2,p_C=>1/3,p_G=>1/8,p_T=>1/24}))
rank jacIpiId

rank jacIpiId_{0,1,2,3,4,5,6,7,8,38} --7
rank jacIpiId_{0,1,2,3,4,5,6,7,8,39} --7
rank jacIpiId_{0,1,2,3,4,5,6,7,8,40} --7
rank jacIpiId_{0,1,2,3,4,5,6,7,8,9,38} --8
for i from 10 to 37 do print(i,rank jacIpiId_{0,1,2,3,4,5,6,7,8,9,i,38}) --8
rank jacIpiId_{0,1,2,3,4,5,6,7,8,9,13,38} --9
rank jacIpiId_{0,1,2,3,4,6,9,13,38} --9
codim ideal {Ipi_0,Ipi_1,Ipi_2,Ipi_3,Ipi_4,Ipi_6,Ipi_9,Ipi_13,Ipi_38} --7
aux=ideal {Ipi_0,Ipi_1,Ipi_2,Ipi_3,Ipi_4,Ipi_6,Ipi_9,Ipi_13,Ipi_38} --7
time minimalPrimes aux
MP=o50
length apply(MP,i->codim i)
aux={(C,C,A),(C,C,G),(T,T,A),(T,T,G),(A,A,A),(G,G,A),(C,C,C),(G,G,G),(T,T,T)}
jacIpiId^(apply(aux,i->position(nonZeroEntries,j->j==i)))

rank jacIpiId_{0,1,2,3,4,5,6,7,8,38,39,40} --7
toList (0..37)
rank jacIpiId_(toList (0..37)) --9
rank jacIpiId_(join({0,1,2,3,4,6},toList (9..37))) --9
rank jacIpiId_(join({0,1,2,3,4,6},toList (9..37))) --9

length toList(9..37)
ind=time subsets(toList(9..37),3);
length ind
rkind=time apply(ind,i->rank jac_(i));
length positions(rkind,i->i==3)
rkindId=time apply(ind,i->rank jacId_(i));
rkindId_(positions(rkindId,i->i!=3))
ind3=set ind - set ind_(positions(rkindId,i->i!=3));

codims=time apply(toList ind3,i->codim ideal apply(i,j->I_j));
length positions(codims,i->i==3)
cod3=(toList ind3)_(positions(codims,i->i==3));
cod3_0
codim ideal{I_24,I_27,I_28}
cod3_1
codim ideal{I_24,I_27,I_29}
netList cod3

CI2=ideal{I2_0,I2_1,I2_2,I2_3,I2_4,I2_6} 
codim CI2--6
codims2=time apply(cod3,i->codim (CI2+ideal apply(i,j->I_j)));
length positions(codims2,i->i==9) 
cod3_(positions(codims2,i->i==9))
CI=CI2+ideal{I_13,I_28,I_37}
codim CI
loadPackage "Depth"
isRegularSequence(gens CI)
netList CI_*

--We can't compute minimal primes in K

CI3=ideal{I_13,I_28,I_37}
ind6=time subsets(toList(0..8),6);
codims6=time apply(toList ind6,i->codim ideal apply(i,j->I_j));
length positions(codims6,i->i==6)
ind6_(positions(codims6,i->i==6))
codims9=time apply(ind6_(positions(codims6,i->i==6)),i->codim (CI3+ideal apply(i,j->I_j)));
length positions(codims9,i->i==9)
codims9
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
mingens I
length gens S
jac=jacobian I;
rank jac 
rank jac_{0,1,2,3,4,5,6,7,8}
det jac_{0,1,2,3,4,5,6,7,8}^{0,1,2,3,4,5,6,7,8}

toString gens S
use S
jacId=sub(jac,{p_(A,A,A)=>p_A, p_(A,C,C)=>0, p_(A,G,G)=>0, p_(A,T,T)=>0, p_(C,A,C)=>0, p_(C,C,A)=>0, p_(C,C,C)=>p_C,
       p_(C,C,G)=>0, p_(C,G,C)=>0, p_(G,A,G)=>0, p_(G,C,C)=>0, p_(G,G,A)=>0, p_(G,G,G)=>p_G, p_(G,T,T)=>0,
       p_(T,A,T)=>0, p_(T,G,T)=>0, p_(T,T,A)=>0, p_(T,T,G)=>0, p_(T,T,T)=>p_T})
rank jacId_{0,1,2,3,4,5,6,7,8} --6
rank jacId                     --19

QBAR=matrix {{p_A+p_C+p_G+p_T}, {(p_C+p_T)/(p_C*p_T)}, {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(p_A+p_G)/(p_A*p_G)}, {(p_C+p_T)/(p_C*p_T)}, {(p_C+p_T)/(p_C*p_T)},
      {(-p_C^2+p_T^2)/(p_C^2*p_T^2)}, {(-1)/(p_C*p_T)}, {(-1)/(p_C*p_T)}, {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(-1)/(p_C*p_T)},
      {(p_A+p_C+p_G+p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)}, {(-p_A^2+p_C^2-2*p_A*p_G-p_G^2+2*p_C*p_T+p_T^2)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C
      *p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2)}, {1/(p_A*p_G)}, {(p_A+p_G)/(p_A*p_G)}, {1/(p_A*p_G)}, {(p_A+p_G)/(p_A*p_G)}, {1/(p_A*p_G)}, {(-p_A^2+p_G^2)/(p_A^2*p_G^2)}}

netList flatten entries QBAR
length unique(flatten entries QBAR)
netList unique(flatten entries QBAR)

jacId=sub(jac,transpose QBAR)
time rank jacId --9
rank jacId_{0,1,2,3,4,5,6,7,8} --6
jacId_{0,1,2,3,4,5,6,7,8}

min9Id=time minors(9,jacId);
--takes more than 324 sec

cols23=toList apply(toList(0..37),i->jacId_i);
length unique cols23


codim ideal{I_0,I_2,I_4,I_5,I_7,I_8}
rank jacId_{0,2,4,5,7,8} --6


min9=time minors(9,jac_{0,1,2,3,4,5,6,7,8});
jac_(toList{9..37})
jac_{38,39,40}
rank jac_{0,1,2,3,4,6,37,38}--8

rank jac_{0,1,2,3,4,5,6,7,8}^{1,2,3,4,5,7,8,9,10}--8

-----------------------------------------------------------------------------
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-------------------------------------------------------------------------------
--Checks with original basis
--uniform distribution
uni=sub(pp,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4});
--b=c=e
EI=sub(pp,flatten toList apply(1..3,i->{c_i=>b_i,e_i=>b_i}));
nonrepeatedEI=time unique(flatten entries EI);
length nonrepeatedEI --44=64-20 
--matches Example pg 1119 in Casanellas-Steel
--but we need 16 more equations!!!
p_T*EI_(position(S,i->i==(A,A,G)),0)==p_G*EI_(position(S,i->i==(A,A,T)),0) 
p_T*EI_(position(S,i->i==(A,A,A)),0)==p_G*EI_(position(S,i->i==(A,A,T)),0) --false
p_T*EI_(position(S,i->i==(A,G,A)),0)==p_G*EI_(position(S,i->i==(A,T,A)),0) 
p_T*EI_(position(S,i->i==(G,A,A)),0)==p_G*EI_(position(S,i->i==(T,A,A)),0) 
p_T*EI_(position(S,i->i==(A,C,G)),0)==p_G*EI_(position(S,i->i==(A,C,T)),0) 
EI_(position(S,i->i==(A,C,G)),0)==EI_(position(S,i->i==(G,A,C)),0) 
p_C*p_T*EI_(position(S,i->i==(A,A,G)),0)==p_A*p_A*EI_(position(S,i->i==(C,G,T)),0)--false


--b=c=e and uniform distr
JC=sub(uni,flatten toList apply(1..3,i->{c_i=>b_i,e_i=>b_i}));
nonrepeatedJC=time unique(flatten entries JC);
length nonrepeatedJC --5=64-59 matches Casanellas-Fernandez-Kedzierska
JC_(position(S,i->i==(A,A,A)),0)==JC_(position(S,i->i==(T,T,T)),0) --true
JC_(position(S,i->i==(A,A,G)),0)==JC_(position(S,i->i==(A,A,T)),0) --true

--c=e and uniform distr
K80=sub(uni,flatten toList apply(1..3,i->{e_i=>c_i}));
nonrepeatedK80=time unique(flatten entries K80);
length nonrepeatedK80 --10=64-54 matches Casanellas-Fernandez-Kedzierska

----------------------------------------------------------------------------
----------------------------------------------------------------------------
-- Checks with original basis rescaled with inverse of pi's appearing in the configuration
qq=mutableMatrix pp;
for i to 63 do (
	qq_(i,0)=1/(p_((S_i)_0)*p_((S_i)_1)*p_((S_i)_2))*qq_(i,0);
	sum apply(OE,k->p_k*M1_(position(OE,l->l==k),position(OE,l->l==(S_i)_0))
	                           *M2_(position(OE,l->l==k),position(OE,l->l==(S_i)_1))
		                   *M3_(position(OE,l->l==k),position(OE,l->l==(S_i)_2)))
) 
qq=matrix qq;
nonrepeatedTN=time unique(flatten entries qq);
length nonrepeatedTN --34 --> 64-34=30 equalities

qq_(position(S,i->i==(A,C,T)),0)==qq_(position(S,i->i==(A,T,C)),0) --true
pp_(position(S,i->i==(C,C,T)),0)==pp_(position(S,i->i==(T,C,C)),0) --false
qq_(position(S,i->i==(A,C,T)),0)==qq_(position(S,i->i==(G,C,T)),0) --true
qq_(position(S,i->i==(G,C,T)),0)==qq_(position(S,i->i==(G,T,C)),0) --true

qq_(position(S,i->i==(A,C,T)),0)==qq_(position(S,i->i==(G,C,T)),0) --true

qq_(position(S,i->i==(C,G,T)),0)==qq_(position(S,i->i==(T,G,C)),0)--true
qq_(position(S,i->i==(C,A,T)),0)==qq_(position(S,i->i==(T,A,C)),0)--true
qq_(position(S,i->i==(C,G,A)),0)==qq_(position(S,i->i==(C,A,G)),0)--true
qq_(position(S,i->i==(C,C,T)),0)==qq_(position(S,i->i==(T,C,C)),0)--false
qq_(position(S,i->i==(C,T,T)),0)==qq_(position(S,i->i==(T,T,C)),0)--false
qq_(position(S,i->i==(C,C,T)),0)==qq_(position(S,i->i==(T,T,C)),0)--false
p_T*qq_(position(S,i->i==(C,C,T)),0)==p_C*qq_(position(S,i->i==(T,T,C)),0)--false

--b=c=e
EI=sub(qq,flatten toList apply(1..3,i->{c_i=>b_i,e_i=>b_i}));
nonrepeatedEI=time unique(flatten entries EI);
length nonrepeatedEI --17 --> 64-17=47 equalities

--c=e
HKY=sub(qq,flatten toList apply(1..3,i->{e_i=>c_i}));
nonrepeatedHKY=time unique(flatten entries HKY);
length nonrepeatedHKY --34

HKY_(position(S,i->i==(A,C,T)),0)==HKY_(position(S,i->i==(A,T,C)),0) 
HKY_(position(S,i->i==(A,C,T)),0)==HKY_(position(S,i->i==(G,C,T)),0) 


