restart

--Ring definition for QUARTETS on p1..p4 and EIGENVALUES l1..l4
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,4)]

--Retrieve tensor on our desired basis
qbar=value get "4leaves_tensor_id_FinalBasis.txt";
pbar=value get "4leaves_tensor_FinalBasis.txt";

--Tensor coordinate indexes
S=sort elements (set {1,2,3,4})^**4/splice/splice;
--Column and row index in flattening matrix
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}

--Build flattening matrix
flattP=mutableMatrix(R,16,16);
for i to 15 do (for j to 15 do
    flattP_(i,j)=pbar_(position(S,k->k==((s_i)_0,(s_i)_1,(s_j)_0,(s_j)_1)),0);
    );
flattP=matrix flattP;

-------------------------------------------------------------
--Invariants of quartets: setup
-------------------------------------------------------------
var=toList apply(S,i->(symbol x)_i)
Rx=K[var]
f=map(R,Rx,flatten entries pbar);
--Example
f(x_(1,1,1,1))
f(x_(1,1,1,2))

-------------------------------------------------------------
--Invariants of quartets arising from invariants of tripods
-------------------------------------------------------------

orderedVars={(1,1,1),(2,2,2),(3,3,3),(4,4,4),(1,2,2),(2,1,2),(2,2,1),(1,3,3),(3,1,3),(3,3,1),(1,4,4),(4,1,4),(4,4,1),(2,3,3),(3,2,3),(3,3,2),(2,4,4),(4,2,4),(4,4,2)}
--netList orderedVars
--length orderedVars --19
varp=toList apply(orderedVars,i->(symbol p)_i);
Rtripod=K[varp]
I=value get "3leavesVanishingIdeal_FinalBasis.txt"

CI=ideal{I_0,I_2,I_6,I_1,I_3,I_7,I_27,I_23,I_13}
codim CI
netList CI_*

s4=apply(orderedVars,i->append(i,1))
s1=apply(orderedVars,i->sequence 1|i);

g4=map(Rx,Rtripod,toList apply(s4,i->x_i));
CI4=g4(CI);
netList CI4_*
apply(flatten entries gens CI4,i->f(i))

g1=map(Rx,Rtripod,toList apply(s1,i->x_i));
CI1=g1(CI);
netList CI1_*
apply(flatten entries gens CI1,i->f(i))

CItripod=trim(CI1+CI4); 
time codim CItripod --14 => not a complete intersection
-- used 174.3 seconds
time codim CI1 --9
 -- used 0.554924 seconds
time codim CI4 --9
  -- used 0.529536 seconds
betti CItripod --18

----------------------------------------
--Egde invariants of quartets for 12|34
----------------------------------------

--Rank-1 invariants:
B1T=flattP_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,2))}
rank B1T

R1T=B1T^(positions(apply(0..15,i->B1T^{i}),j->j!=0))
s_(positions(apply(0..15,i->B1T^{i}),j->j!=0))

          14        j
14     (14,14)   (14,j)
i       (i,14)    (i,j)

rk1T=flatten toList apply({(4,1),(2,4),(4,2)},
    j->apply({(4,1),(2,4),(4,2),(4,4)},i->x_(1,4,1,4)*x_(i|j)-x_((1,4)|j)*x_(i|(1,4))))
length rk1T --12

B1C=flattP_{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(2,3)),position(s,i->i==(3,2))}
rank B1C

R1C=B1C^(positions(apply(0..15,i->B1C^{i}),j->j!=0))
s_(positions(apply(0..15,i->B1C^{i}),j->j!=0))

         13         j
13    (13,13)    (13,j)
i      (i,13)     (i,j)

rk1C=flatten toList apply({(3,1),(2,3),(3,2)},j->apply({(3,1),(2,3),(3,2),(3,3)},i->x_(1,3,1,3)*x_(i|j)-x_((1,3)|j)*x_(i|(1,3))))
length rk1C --12

--Double-check
apply(rk1C,i->f(i))
apply(rk1T,i->f(i))
f(x_(1,3,3,1)*x_(3,1,1,3))==f(x_(1,3,1,3)*x_(3,1,3,1))


B1G=flattP_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}
rank B1G

R1G=B1G^(positions(apply(0..15,i->B1G^{i}),j->j!=0))
s_(positions(apply(0..15,i->B1G^{i}),j->j!=0))

R1G=flattp_{position(s,i->i==(A,G)),position(s,i->i==(G,A))}^{position(s,i->i==(T,T)),position(s,i->i==(G,G)),position(s,i->i==(A,G)),position(s,i->i==(G,A)),position(s,i->i==(C,C))}

         12         21
12    (12,12)     (12,21)
i      (i,12)      (i,21)

rk1G=toList apply({(2,1),(3,3),(2,2),(4,4)},i->x_(1,2,1,2)*x_(i|(2,1))-x_(1,2,2,1)*x_(i|(1,2)))
length rk1G --4
--Double-check
apply(rk1G,i->f(i))

rk1=rk1C|rk1G|rk1T
length rk1

-- Rank-2 invariants

B2G=flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}
rank B2G

R2G=B2G^(positions(apply(0..15,i->B2G^{i}),j->j!=0))
s_(positions(apply(0..15,i->B2G^{i}),j->j!=0))


        11     12     22
   
11     1111   1112   1122

12     1211   1212   1222

i      i11     i12   i22


rk2=toList apply({(2,1),(3,3),(2,2),(4,4)},i->det matrix{{x_(1,1,1,1),x_(1,1,1,2),x_(1,1,2,2)},
	{x_(1,2,1,1),x_(1,2,1,2),x_(1,2,2,2)},{x_(i|(1,1)),x_(i|(1,2)),x_(i|(2,2))}})

apply(rk2,i->f(i))
length rk2 --4

-- Rank-3 invariants

B3T=flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}
rank B3T

R3T=B3T^(positions(apply(0..15,i->B3T^{i}),j->j!=0))
s_(positions(apply(0..15,i->B3T^{i}),j->j!=0))
toString oo
{(1,1), (1,4), (4,1), (2,4), (4,2), (4,4), (2,2), (1,2), (2,1), (3,3)}

        11     12     14    44
   
11     1111   1112   1114  1144

12     1211   1212   1214  1244

14     1411   1412   1414  1444

i      i11    i12    i14   i44

rk3T=toList apply({(4,1),(2,4),(4,2),(4,4),(2,2),(2,1),(3,3)},
    i->det matrix{{x_(1,1,1,1),x_(1,1,1,2),x_(1,1,1,4),x_(1,1,4,4)},
	          {x_(1,2,1,1),x_(1,2,1,2),x_(1,2,1,4),x_(1,2,4,4)},
	          {x_(1,4,1,1),x_(1,4,1,2),x_(1,4,1,4),x_(1,4,4,4)},
	          {x_(i|(1,1)),x_(i|(1,2)),x_(i|(1,4)),x_(i|(4,4))}})

apply(rk3T,i->f(i))
length rk3T --7

B3C=flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,3))}
rank B3C

R3C=B3C^(positions(apply(0..15,i->B3C^{i}),j->j!=0))
s_(positions(apply(0..15,i->B3C^{i}),j->j!=0))
toString oo
{(1,1), (4,4), (2,2), (1,2), (2,1), (3,3), (1,3), (3,1), (2,3), (3,2)}

        11     12     13    33
   
11     1111   1112   1113  1133

12     1211   1212   1213  1233

13     1311   1312   1313  1333

i      i11    i12    i13   i33

rk3C=toList apply({(4,4),(2,2),(2,1),(3,3),(3,1),(2,3),(3,2)},
       i->det matrix{{x_(1,1,1,1),x_(1,1,1,2),x_(1,1,1,3),x_(1,1,3,3)},
	             {x_(1,2,1,1),x_(1,2,1,2),x_(1,2,1,3),x_(1,2,3,3)},
	             {x_(1,3,1,1),x_(1,3,1,2),x_(1,3,1,3),x_(1,3,3,3)},
	             {x_(i|(1,1)),x_(i|(1,2)),x_(i|(1,3)),x_(i|(3,3))}})

apply(rk3C,i->f(i))
length rk3C

rk3=rk3C|rk3T
length rk3 --14

--Rank-0 invariants

rk0={x_(3,4,3,4),x_(3,4,4,3),x_(4,3,3,4),x_(4,3,4,3)}
apply(rk0,i->f(i))

edge=rk0|rk1|rk2|rk3
CIedge=trim ideal edge;
CI=CIedge+CItripod;
betti CI --68
numgens CI  --18 
numgens CIedge --50

jac=jacobian CI;
time rank jac --68 
--here we don't really know if all 68-minors could be 0

--no evolution point
p0=flatten entries sub(sub(sub(qbar,R),apply(gens R,i->i=>1)),Rx);
jacId=sub(jac,matrix{p0});
time rank jacId --68

mins=time minors(68,jacId);
mins
