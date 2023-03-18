-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
--------------------- HOMOGENEOUS PARAMETRIZATION FOR TN93 --------------------------------------------
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
Rl=K[l_(5,A),l_(5,C),l_(5,G),l_(5,T)]
--basis
--Hl=transpose(matrix{{4,4,4,4},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
Hl=transpose(matrix{{1,1,1,1},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});

--substitution matrix in the eigenvalues
Ml=Hl*diagonalMatrix(Rl,4,4,{l_(5,A),l_(5,C),l_(5,G),l_(5,T)})*inverse(Hl)
inverse(Hl)*Ml*Hl
--Flattening 12|34 for identity at the leaves
OE={A,C,G,T}
flattql=mutableMatrix(Rl,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattql_(i,j)=p_(OE_(i%4))*Ml_(i%4,j%4));      
    );
flattql=matrix(flattql);
--Change of basis
flattQl=time (transpose(Hl)**transpose(Hl))*flattql*(Hl**Hl);
--Quasi-block form
blockQl=flattQl_{0,3,12,14,11,15,10,2,8,5,1,4,9,6,13,7}^{0,3,12,14,11,15,10,2,8,5,1,4,9,6,13,7};
blockQl_{0,1,5,9}

rank blockQl_{0,1,5,9}
blockId=sub(blockQl,{l_(5,A)=>1,l_(5,C)=>1,l_(5,G)=>1,l_(5,T)=>1});
rank blockId
rank blockId_{0,1,5,9}
rank blockId_{0,1,5,9}^{0,5,6,9}
blockId_{5,6,9}^{5,6,9}

-------------------------------
--Examples for a given pi
-------------------------------
Rex=K[lA,lC,lG,lT];
f=map(Rex,Rl,gens Rex);
blockEx=f(blockQl);
sub(blockEx,{p_A=>1/2,p_C=>1/3,p_G=>1/8,p_T=>1/24})
sub(blockEx,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4})

------------------------------
-- Checks
------------------------------
--Check homogeneity
isHomogeneous(ideal flatten entries blockQl) --true
betti(ideal flatten entries blockQl) --linear in lambda (80 non-zero entries)
--Check equality with parametrization in overleaf
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

--It's enough to assume that p_A+p_C+p_G+p_T=1 to observe equality:
sub(blockQl-pTN92,{p_A=>1-p_C-p_G-p_T}) --0
---------------------------------------------------------------------------
----------------------------------------------------------------------------

--Ring for general quartet
Rgeneral=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T),l_(4,A),l_(4,C),l_(4,G),l_(4,T),l_(5,A),l_(5,C),l_(5,G),l_(5,T)]
--Flattening 12|34 for general quartet
s={(A,A),(A,T),(T,A),(T,G),(G,T),(T,T),(G,G),(A,G),(G,A),(C,C),(A,C),(C,A),(C,G),(G,C),(T,C),(C,T)}
Q=sub(blockQl,Rgeneral);
q=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*Q_(i,j)));
--Example
blockQl_{0,1,2,3,4}^{0,1,2,3,4}
q_{0,1,2,3,4}^{0,1,2,3,4}
sub(q,{p_A=>1/2,p_C=>1/3,p_G=>1/8,p_T=>1/24})
rank q_{1,2,3,4}^{1,2,3,4}

--Check homogeneity
isHomogeneous(ideal flatten entries q) --true
betti(ideal flatten entries q) --quintic in lambda (80 non-zero entries)


varp=flatten toList apply(0..15,i->toList apply(0..15,j->(symbol p)_(s_i,s_j)));
S=QQ[varp];
length gens S

P=transpose genericMatrix(S,p_((A,A),(A,A)),16,16)

degp=apply(gens S,i->5)
degl=apply(gens Rgeneral,i->1)
deg=join(degp,degl)

RS=K[gens S,gens Rgeneral,Degrees=>deg]
P=sub(P,RS);
q=sub(q,RS);
I=trim ideal flatten entries (P-q);
betti I
isHomogeneous I --true

toExternalString I
toString I

--Ideal definition for TN93
J=time eliminate(apply(gens Rgeneral,i->sub(i,RS)),I);

