-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
------------ Parametrization of TN92 in parameters b,c,e ----------------------
------------------ with identity at the leaves --------------------------------
-------------------------------------------------------------------------------

restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
R=K[b_5,c_5,e_5];
i=5;
m=matrix{{0,b_i,c_i,b_i},{b_i,0,b_i,e_i},{c_i,b_i,0,b_i},{b_i,e_i,b_i,0}}

OE={A,C,G,T}
mm=mutableMatrix {toList apply(OE,i->p_(OE_((position(OE,j->j==i))%4))*(matrix m)_{position(OE,j->j==i)})}
for i to 3 do mm_(i,i)=1-sum flatten entries (matrix mm)^{i}
M=matrix mm
for i to 3 do print sum flatten entries M^{0}
for i to 3 do print sum flatten entries M_{0}

--Basis that provides many 0's and has signs correponding to Hadamart matrix
H=transpose(matrix{{4,4,4,4},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
 
Hinv=inverse H
D=Hinv*M*H

--parameters b,c,e can be expressed in terms of eigenvalues

1/(p_A+p_C+p_G+p_T)*(D_(0,0)-D_(2,2))

1/((p_A+p_G)*(p_A+p_C+p_G+p_T))*((p_A+p_G)*D_(0,0)+(p_C+p_T)*D_(2,2)-(p_A+p_C+p_G+p_T)*D_(3,3))

1/((p_C+p_T)*(p_A+p_C+p_G+p_T))*((p_C+p_T)*D_(0,0)+(p_A+p_G)*D_(2,2)-(p_A+p_C+p_G+p_T)*D_(1,1))

--Flattening
flattq=mutableMatrix(R,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=p_(OE_(i%4))*M_(i%4,j%4));      
    )
flattq=matrix(flattq)

--Change of basis
flattQ=time (transpose(H)**transpose(H))*flattq*(H**H);

--Quasi-block form
blockQ=flattQ_{0,3,12,14,11,15,10,2,8,5,1,4,9,6,13,7}^{0,3,12,14,11,15,10,2,8,5,1,4,9,6,13,7};

B1=blockQ^{1,2,3,4}_{1,2,3,4}
rank B1 
B2=blockQ^{6,7,8,9}_{6,7,8,9}
rank B2 
B3=blockQ^{10,11,12,13}_{10,11,12,13}
rank B3 

--Substitute by certain values
UD=sub(blockQ,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4});
UD^{1,2,3,4}_{1,2,3,4}
UD

UDJC=sub(UD,{c_5=>b_5,e_5=>b_5});
UDJC^{1,2,3,4}_{1,2,3,4}
UDJC

rank UDJC --4

qtest2=sub(blockQ,{p_A=>1/3,p_C=>1/3,p_G=>1/6,p_T=>1/6,c_5=>b_5,e_5=>b_5});
rank qtest2

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
------------------- Parametrization of TN92 in the eigenvalues ----------------
---------------------- with identity at the leaves ----------------------------
-------------------------------------------------------------------------------
Rl=K[l_(5,A),l_(5,C),l_(5,G),l_(5,T)]
f=map(Rl,R,{1/(p_A+p_C+p_G+p_T)*(l_(5,A)-l_(5,G)),
1/((p_A+p_G)*(p_A+p_C+p_G+p_T))*((p_A+p_G)*l_(5,A)+(p_C+p_T)*l_(5,G)-(p_A+p_C+p_G+p_T)*l_(5,T)),
1/((p_C+p_T)*(p_A+p_C+p_G+p_T))*((p_C+p_T)*l_(5,A)+(p_A+p_G)*l_(5,G)-(p_A+p_C+p_G+p_T)*l_(5,C))})

--Check values of b,c,e depending on pi
L={1/(p_A+p_C+p_G+p_T)*(l_(5,A)-l_(5,G)),
1/((p_A+p_G)*(p_A+p_C+p_G+p_T))*((p_A+p_G)*l_(5,A)+(p_C+p_T)*l_(5,G)-(p_A+p_C+p_G+p_T)*l_(5,T)),
1/((p_C+p_T)*(p_A+p_C+p_G+p_T))*((p_C+p_T)*l_(5,A)+(p_A+p_G)*l_(5,G)-(p_A+p_C+p_G+p_T)*l_(5,C))}
apply(L,i->sub(i,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4}))
apply(L,i->sub(i,{p_A=>1/3,p_C=>1/3,p_G=>1/6,p_T=>1/6}))
apply(L,i->sub(i,{p_A=>1/2,p_C=>1/3,p_G=>1/8,p_T=>1/24}))
--uniform dist + 3 eigenvalues equal
apply(L,i->sub(i,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4,l_(5,A)=>1,l_(5,G)=>l_(5,C),l_(5,T)=>l_(5,C)}))
--> b=c=e

--Recompute matrices in the lambdas
Ml=f(M);
Hl=sub(H,Rl);
Dl=(inverse Hl)*Ml*Hl

--Flattening of f(M)
flattql=mutableMatrix(Rl,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattql_(i,j)=p_(OE_(i%4))*Ml_(i%4,j%4));      
    );
flattql=matrix(flattql);

--Change of basis
flattQl=time (transpose(Hl)**transpose(Hl))*flattql*(Hl**Hl);

--Quasi-block form
blockQl=flattQl_{0,3,12,14,11,15,10,2,8,5,1,4,9,6,13,7}^{0,3,12,14,11,15,10,2,8,5,1,4,9,6,13,7};

B1l=blockQl^{1,2,3,4}_{1,2,3,4}
rank B1l 
B2l=blockQl^{6,7,8,9}_{6,7,8,9}
rank B2l 
B3l=blockQl^{10,11,12,13}_{10,11,12,13}
rank B3l 

--Study degrees and monomiality/binomiality of the entries
length (flatten entries blockQl)
param=select(unique(flatten entries blockQl),i->i!=0)
length param
netList param
netList apply(param,i->sub(i,{p_A=>1-p_C-p_G-p_T})) --> doesn't make it better
deg0=select(param,i->degree(i)=={0})
deg1=select(param,i->degree(i)=={1})
length deg0,length deg1

netList apply(deg1,i->support i)
position(deg1,i->select(terms i,i->degree(i)=={0})=={})
netList terms deg1_8
list0=flatten toList apply(deg1,i->select(terms i,i->degree(i)=={0}))
listA=flatten toList apply(deg1,i->select(terms i,i->support(i)=={l_(5,A)}))
length list0
length listA
netList list0
netList drop(listA,{8,8})

netList apply(list0,drop(listA,{8,8}),(i,j)->i*l_(5,A)+j)
netList apply(list0,drop(deg1,{8,8}),(i,j)->length terms(j-i+i*l_(5,A)))
--If we assume that l_A=1, then most degree 1 entries are monomial in the lambdas. Only 6 aren't
--5 in this list plus deg1_8

--Comparison with previous parametrization
-- Summary: when assuming l_A=1 and sum p_i=1, it's exactly the same parametrization
pTN92=matrix {{256*l_(5,A), 0, 0, 0, 0, ((16*p_A+16*p_G)/(p_A*p_G))*l_(5,A), (16/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(5,A), 0, 0, ((16*p_C+16*p_T)/(p_C*p_T))*l_(5,A), 0, 0, 0, 0, 0, 0},
{0, ((16*p_A+16*p_G)/(p_A*p_G))*l_(5,T), ((16*p_A+16*p_G)/(p_A*p_G))*l_(5,T), (4/(p_A*p_G))*l_(5,T), (4/(p_A*p_G))*l_(5,T), ((-4*p_A^2+4*p_G^2)/(p_A^2*p_G^2))*l_(5,T), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
{0, ((16*p_A+16*p_G)/(p_A*p_G))*l_(5,T), ((16*p_A+16*p_G)/(p_A*p_G))*l_(5,T), (4/(p_A*p_G))*l_(5,T), (4/(p_A*p_G))*l_(5,T), ((-4*p_A^2+4*p_G^2)/(p_A^2*p_G^2))*l_(5,T),0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
{0, (4/(p_A*p_G))*l_(5,T), (4/(p_A*p_G))*l_(5,T), (1/(p_A^2*p_G+p_A*p_G^2))*l_(5,T),(1/(p_A^2*p_G+p_A*p_G^2))*l_(5,T),((-p_A+p_G)/(p_A^2*p_G^2))*l_(5,T), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
{0, (4/(p_A*p_G))*l_(5,T), (4/(p_A*p_G))*l_(5,T), (1/(p_A^2*p_G+p_A*p_G^2))*l_(5,T),(1/(p_A^2*p_G+p_A*p_G^2))*l_(5,T), ((-p_A+p_G)/(p_A^2*p_G^2))*l_(5,T), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
{((16*p_A+16*p_G)/(p_A*p_G))*l_(5,A), ((-4*p_A^2+4*p_G^2)/(p_A^2*p_G^2))*l_(5,T),
      ((-4*p_A^2+4*p_G^2)/(p_A^2*p_G^2))*l_(5,T), ((-p_A+p_G)/(p_A^2*p_G^2))*l_(5,T), ((-p_A+p_G)/(p_A^2*p_G^2))*l_(5,T),
      ((-p_A*p_C-p_C*p_G-p_A*p_T-p_G*p_T+p_A+p_G)/(p_A^2*p_G^2))*l_(5,A)+((p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)/(p_A^2*p_G^2))*l_(5,G)+((p_A^3-p_A^2*p_G-p_A*p_G^2+p_G^3)/(p_A^3*p_G^3))*l_(5,T),
      (1/(p_A*p_C*p_G+p_A*p_G*p_T))*l_(5,A)+((-p_A+p_C-p_G+p_T)/(p_A^2*p_C*p_G+p_A*p_C*p_G^2+p_A^2*p_G*p_T+p_A*p_G^2*p_T))*l_(5,G), (4/(p_A*p_G))*l_(5,G), (4/(p_A*p_G))*l_(5,G),
      ((p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)/(p_A*p_C*p_G*p_T))*l_(5,A)+((-p_A*p_C-p_C*p_G-p_A*p_T-p_G*p_T)/(p_A*p_C*p_G*p_T))*l_(5,G), 0, 0, 0, 0, 0, 0},
      {(16/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(5,A), 0, 0, 0, 0,
      (1/(p_A*p_C*p_G+p_A*p_G*p_T))*l_(5,A)+((-p_A+p_C-p_G+p_T)/(p_A^2*p_C*p_G+p_A*p_C*p_G^2+p_A^2*p_G*p_T+p_A*p_G^2*p_T))*l_(5,G),
      ((p_A^3+p_C^3+3*p_A^2*p_G+3*p_A*p_G^2+p_G^3+3*p_C^2*p_T+3*p_C*p_T^2+p_T^3-4*p_C^2-8*p_C*p_T-4*p_T^2+4*p_C+4*p_T-1)/(p_A^3*p_C^3+3*p_A^2*p_C^3*p_G+3*p_A*p_C^3*p_G^2+p_C^3*p_G^3+3*p_A^3
      *p_C^2*p_T+9*p_A^2*p_C^2*p_G*p_T+9*p_A*p_C^2*p_G^2*p_T+3*p_C^2*p_G^3*p_T+3*p_A^3*p_C*p_T^2+9*p_A^2*p_C*p_G*p_T^2+9*p_A*p_C*p_G^2*p_T^2+3*p_C*p_G^3*p_T^2+p_A^3*p_T^3+3*p_A^2*p_G*p_T^3+
      3*p_A*p_G^2*p_T^3+p_G^3*p_T^3))*l_(5,A)+((4*p_C^2+8*p_C*p_T+4*p_T^2-4*p_C-4*p_T+1)/(p_A^3*p_C^3+3*p_A^2*p_C^3*p_G+3*p_A*p_C^3*p_G^2+p_C^3*p_G^3+3*p_A^3*p_C^2*p_T+9*p_A^2*p_C^2*p_G*p_T
      +9*p_A*p_C^2*p_G^2*p_T+3*p_C^2*p_G^3*p_T+3*p_A^3*p_C*p_T^2+9*p_A^2*p_C*p_G*p_T^2+9*p_A*p_C*p_G^2*p_T^2+3*p_C*p_G^3*p_T^2+p_A^3*p_T^3+3*p_A^2*p_G*p_T^3+3*p_A*p_G^2*p_T^3+p_G^3*p_T^3))*
      l_(5,G), ((-4*p_A+4*p_C-4*p_G+4*p_T)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(5,G),
      ((-4*p_A+4*p_C-4*p_G+4*p_T)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(5,G),
      (1/(p_A*p_C*p_T+p_C*p_G*p_T))*l_(5,A)+((p_A-p_C+p_G-p_T)/(p_A*p_C^2*p_T+p_C^2*p_G*p_T+p_A*p_C*p_T^2+p_C*p_G*p_T^2))*l_(5,G), 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, (4/(p_A*p_G))*l_(5,G),
      ((-4*p_A+4*p_C-4*p_G+4*p_T)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(5,G),
      ((16*p_A+16*p_C+16*p_G+16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(5,G), ((16*p_A+16*p_C+16*p_G+16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(5,G), ((-4)/(p_C*p_T))*l_(5,G), 0, 0,
      0, 0, 0, 0}, {0, 0, 0, 0, 0, (4/(p_A*p_G))*l_(5,G), ((-4*p_A+4*p_C-4*p_G+4*p_T)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+
      2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(5,G), ((16*p_A+16*p_C+16*p_G+16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(5,G),
      ((16*p_A+16*p_C+16*p_G+16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(5,G), ((-4)/(p_C*p_T))*l_(5,G), 0, 0, 0, 0, 0, 0}, {((16*p_C+16*p_T)/(p_C*p_T))*l_(5,A), 0, 0, 0, 0,
      ((p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)/(p_A*p_C*p_G*p_T))*l_(5,A)+((-p_A*p_C-p_C*p_G-p_A*p_T-p_G*p_T)/(p_A*p_C*p_G*p_T))*l_(5,G),
      (1/(p_A*p_C*p_T+p_C*p_G*p_T))*l_(5,A)+((p_A-p_C+p_G-p_T)/(p_A*p_C^2*p_T+p_C^2*p_G*p_T+p_A*p_C*p_T^2+p_C*p_G*p_T^2))*l_(5,G), ((-4)/(p_C*p_T))*l_(5,G), ((-4)/(p_C*p_T))*l_(5,G),
      ((-p_A*p_C-p_C*p_G-p_A*p_T-p_G*p_T+p_C+p_T)/(p_C^2*p_T^2))*l_(5,A)+((p_C^3-p_C^2*p_T-p_C*p_T^2+p_T^3)/(p_C^3*p_T^3))*l_(5,C)+((p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)/(p_C^2*p_T^2))*l_(5,G),
      ((-4*p_C^2+4*p_T^2)/(p_C^2*p_T^2))*l_(5,C), ((-4*p_C^2+4*p_T^2)/(p_C^2*p_T^2))*l_(5,C), ((p_C-p_T)/(p_C^2*p_T^2))*l_(5,C), ((p_C-p_T)/(p_C^2*p_T^2))*l_(5,C), 0, 0}, {0, 0, 0, 0, 0, 0,
      0, 0, 0, ((-4*p_C^2+4*p_T^2)/(p_C^2*p_T^2))*l_(5,C), ((16*p_C+16*p_T)/(p_C*p_T))*l_(5,C), ((16*p_C+16*p_T)/(p_C*p_T))*l_(5,C), ((-4)/(p_C*p_T))*l_(5,C), ((-4)/(p_C*p_T))*l_(5,C), 0,
      0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, ((-4*p_C^2+4*p_T^2)/(p_C^2*p_T^2))*l_(5,C), ((16*p_C+16*p_T)/(p_C*p_T))*l_(5,C), ((16*p_C+16*p_T)/(p_C*p_T))*l_(5,C), ((-4)/(p_C*p_T))*l_(5,C),
      ((-4)/(p_C*p_T))*l_(5,C), 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, ((p_C-p_T)/(p_C^2*p_T^2))*l_(5,C), ((-4)/(p_C*p_T))*l_(5,C), ((-4)/(p_C*p_T))*l_(5,C), (1/(p_C^2*p_T+p_C*p_T^2))*l_(5,C),
      (1/(p_C^2*p_T+p_C*p_T^2))*l_(5,C), 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, ((p_C-p_T)/(p_C^2*p_T^2))*l_(5,C), ((-4)/(p_C*p_T))*l_(5,C), ((-4)/(p_C*p_T))*l_(5,C),
      (1/(p_C^2*p_T+p_C*p_T^2))*l_(5,C), (1/(p_C^2*p_T+p_C*p_T^2))*l_(5,C), 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

aux=blockQl-pTN92;
aux_{0,1,2,3}^{0,1,2,3}
auxA1=sub(aux,{l_(5,A)=>1});
auxA1_{0,1,2,3}^{0,1,2,3}
netList select(unique flatten entries auxA1,i->i!=0)
sub(auxA1,{p_A=>1-p_C-p_G-p_T}) --enough to get zero
sub(aux,{p_A=>1-p_C-p_G-p_T}) --not enough to get zero

--Tests:

--Uniform distribution
test1=sub(blockQl,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4})
--Uniform distribution + 3 equal eigenvalues
test2=sub(blockQl,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4,l_(5,A)=>1,l_(5,G)=>l_(5,C),l_(5,T)=>l_(5,C)})
rank test2
test3=sub(blockQl,{p_A=>1/2,p_C=>1/3,p_G=>1/8,p_T=>1/24})
toString test3
matrix {{256, 0, 0, 0, 0, 160, 1024/15, 0, 0, 432, 0, 0, 0, 0, 0, 0}, {0, -160*l_(5,A)+160*l_(5,T)+160, -160*l_(5,A)+160*l_(5,T)+160, -64*l_(5,A)+64*l_(5,T)+64,
       -64*l_(5,A)+64*l_(5,T)+64, 240*l_(5,A)-240*l_(5,T)-240, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -160*l_(5,A)+160*l_(5,T)+160, -160*l_(5,A)+160*l_(5,T)+160, -64*l_(5,A)+64*l_(5,T)+64,
       -64*l_(5,A)+64*l_(5,T)+64, 240*l_(5,A)-240*l_(5,T)-240, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -64*l_(5,A)+64*l_(5,T)+64, -64*l_(5,A)+64*l_(5,T)+64,
       ((-128)/5)*l_(5,A)+(128/5)*l_(5,T)+128/5, ((-128)/5)*l_(5,A)+(128/5)*l_(5,T)+128/5, 96*l_(5,A)-96*l_(5,T)-96, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -64*l_(5,A)+64*l_(5,T)+64,
       -64*l_(5,A)+64*l_(5,T)+64, ((-128)/5)*l_(5,A)+(128/5)*l_(5,T)+128/5, ((-128)/5)*l_(5,A)+(128/5)*l_(5,T)+128/5, 96*l_(5,A)-96*l_(5,T)-96, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {160,
       240*l_(5,A)-240*l_(5,T)-240, 240*l_(5,A)-240*l_(5,T)-240, 96*l_(5,A)-96*l_(5,T)-96, 96*l_(5,A)-96*l_(5,T)-96, -420*l_(5,A)+60*l_(5,G)+360*l_(5,T)+520,
       (256/15)*l_(5,A)+((-256)/15)*l_(5,G)+128/5, -64*l_(5,A)+64*l_(5,G)+64, -64*l_(5,A)+64*l_(5,G)+64, 270*l_(5,A)-270*l_(5,G), 0, 0, 0, 0, 0, 0}, {1024/15, 0, 0, 0, 0,
       (256/15)*l_(5,A)+((-256)/15)*l_(5,G)+128/5, ((-16384)/3375)*l_(5,A)+(16384/3375)*l_(5,G)+77824/3375, (4096/225)*l_(5,A)+((-4096)/225)*l_(5,G)+(-4096)/225,
       (4096/225)*l_(5,A)+((-4096)/225)*l_(5,G)+(-4096)/225, ((-384)/5)*l_(5,A)+(384/5)*l_(5,G)+192, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -64*l_(5,A)+64*l_(5,G)+64,
       (4096/225)*l_(5,A)+((-4096)/225)*l_(5,G)+(-4096)/225, ((-1024)/15)*l_(5,A)+(1024/15)*l_(5,G)+1024/15, ((-1024)/15)*l_(5,A)+(1024/15)*l_(5,G)+1024/15, 288*l_(5,A)-288*l_(5,G)-288, 0,
       0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -64*l_(5,A)+64*l_(5,G)+64, (4096/225)*l_(5,A)+((-4096)/225)*l_(5,G)+(-4096)/225, ((-1024)/15)*l_(5,A)+(1024/15)*l_(5,G)+1024/15,
       ((-1024)/15)*l_(5,A)+(1024/15)*l_(5,G)+1024/15, 288*l_(5,A)-288*l_(5,G)-288, 0, 0, 0, 0, 0, 0}, {432, 0, 0, 0, 0, 270*l_(5,A)-270*l_(5,G), ((-384)/5)*l_(5,A)+(384/5)*l_(5,G)+192,
       288*l_(5,A)-288*l_(5,G)-288, 288*l_(5,A)-288*l_(5,G)-288, -13122*l_(5,A)+11907*l_(5,C)+1215*l_(5,G)+13851, 2268*l_(5,A)-2268*l_(5,C)-2268, 2268*l_(5,A)-2268*l_(5,C)-2268,
       -1512*l_(5,A)+1512*l_(5,C)+1512, -1512*l_(5,A)+1512*l_(5,C)+1512, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 2268*l_(5,A)-2268*l_(5,C)-2268, -432*l_(5,A)+432*l_(5,C)+432,
       -432*l_(5,A)+432*l_(5,C)+432, 288*l_(5,A)-288*l_(5,C)-288, 288*l_(5,A)-288*l_(5,C)-288, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 2268*l_(5,A)-2268*l_(5,C)-2268,
       -432*l_(5,A)+432*l_(5,C)+432, -432*l_(5,A)+432*l_(5,C)+432, 288*l_(5,A)-288*l_(5,C)-288, 288*l_(5,A)-288*l_(5,C)-288, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0,
       -1512*l_(5,A)+1512*l_(5,C)+1512, 288*l_(5,A)-288*l_(5,C)-288, 288*l_(5,A)-288*l_(5,C)-288, -192*l_(5,A)+192*l_(5,C)+192, -192*l_(5,A)+192*l_(5,C)+192, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0,
       0, -1512*l_(5,A)+1512*l_(5,C)+1512, 288*l_(5,A)-288*l_(5,C)-288, 288*l_(5,A)-288*l_(5,C)-288, -192*l_(5,A)+192*l_(5,C)+192, -192*l_(5,A)+192*l_(5,C)+192, 0, 0}, {0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}


-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
------------------- Parametrization of TN93 in the eigenvalues ----------------
------------------------------------ general ----------------------------------
-------------------------------------------------------------------------------
Rgeneral=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T),l_(4,A),l_(4,C),l_(4,G),l_(4,T),l_(5,A),l_(5,C),l_(5,G),l_(5,T)]

s={(A,A),(A,T),(T,A),(T,G),(G,T),(T,T),(G,G),(A,G),(G,A),(C,C),(A,C),(C,A),(C,G),(G,C),(T,C),(C,T)}
Q=sub(blockQl,Rgeneral);
q=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*Q_(i,j)));

varp=flatten toList apply(0..15,i->toList apply(0..15,j->(symbol p)_(s_i,s_j)));
S=QQ[varp];

P=transpose genericMatrix(S,p_((A,A),(A,A)),16,16)

RS=K[gens S,gens Rgeneral]
P=sub(P,RS);
q=sub(q,RS);
I=trim ideal flatten entries (P-q);
betti I
-- 176 linear equations, 7 quartics, 73 quintics
param=unique(flatten entries gens I);
length param
netList param
apply(param,i->degree i)
netList flatten toList apply(param,i->select(terms i,i->degree(i)=={1}))

deg1=param_(positions(param,i->degree i=={1}));
netList deg1
deg4=param_(positions(param,i->degree i=={4}));
netList deg4
deg5=param_(positions(param,i->degree i=={5}));
netList deg5

--Ideal definition for TN93
J=time eliminate(apply(gens Rgeneral,i->sub(i,RS)),I);

--Ideal definition for 3 equal eigenvalues and pi_i=1/4
netList support I
sub(I,{p_A=>1/4,p_C=>1/4}); --works
sub(I,{p_A=>1/4,p_C=>1/4,p_G=>1/4}); --Process M2 exited abnormally with code 1
IJC=sub(I,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4,l_(5,A)=>1,l_(5,G)=>l_(5,C),l_(5,T)=>l_(5,C)});
