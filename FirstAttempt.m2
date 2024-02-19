restart
K=frac(QQ[p_1,p_2,p_3,p_4,t]);
--R=QQ[p_1,p_2,p_3,p_4,t]/(1-(p_1+p_2+p_3+p_4))
--K=frac(R);
--R=QQ[p_1,p_2,p_3,p_4,t]
--I=trim ideal{(p_1+p_2+p_3+p_4)*t,1-(p_1+p_2+p_3+p_4)}

--M=matrix{{1-(p_2+p_3+p_4)*t,p_2*t,p_3*t,p_4*t},{p_1*t,1-(p_1+p_3+p_4)*t,p_3*t,p_4*t},{p_1*t,p_2*t,1-(p_1+p_2+p_4)*t,p_4*t},{p_1*t,p_2*t,p_3*t,1-(p_1+p_2+p_3)*t}}
M=matrix{{1-(1-p_1)*t,p_2*t,p_3*t,p_4*t},{p_1*t,1-(1-p_2)*t,p_3*t,p_4*t},{p_1*t,p_2*t,1-(1-p_3)*t,p_4*t},{p_1*t,p_2*t,p_3*t,1-(1-p_4)*t}}

flattq=mutableMatrix(K,16,16)
flattq 

for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=p_(i%4+1)*M_(i%4,j%4));      
    )

flattq=matrix(flattq)

--Base Angelica
--H=transpose(matrix{{1,1,1,1},{-p_2,p_1,0,0},{-p_3,0,p_1,0},{-p_4,0,0,p_1}});
--Hinv=inverse(H);

--Base Marta
H=transpose(matrix{{1,1,1,1},{1/p_1,1/p_2,-1/p_3,-1/p_4},{1/p_1,-1/p_2,1/p_3,-1/p_4},{1/p_1,-1/p_2,-1/p_3,1/p_4}})
Hinv=inverse(H);
--In this particular case:
Hinv=(p_1+p_2+p_3+p_4)*Hinv
D=Hinv*M*H
H*M*Hinv
(p_1+p_2+p_3+p_4)*H*M*Hinv

Hinv**Hinv

flattQ=time (Hinv**Hinv)*flattq*transpose(Hinv**Hinv);

numerator flattQ_(0,0)
denominator flattQ_(0,0)

apply(flatten entries flattQ,i->denominator i)
flattQR=sub(flattQ,R);

L=flatten entries flattQR;
Lsim=apply(flatten entries flattQR,i->(trim ideal{i,1-(p_1+p_2+p_3+p_4)})_1);

L_17
Lsim_17
aux=trim ideal{L_17,1-(p_1+p_2+p_3+p_4)}
netList aux_*

--------Jukes-Cantor
restart
K=frac(QQ[t])
M=matrix{{1-3*t,t,t,t},{t,1-3*t,t,t},{t,t,1-3*t,t},{t,t,t,1-3*t}}
flattq=mutableMatrix(K,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=(1/4)*M_(i%4,j%4));      
    )
flattq=matrix(flattq)

H=sub(matrix{{1,1,1,1},{1,1,-1,-1},{1,-1,1,-1},{1,-1,-1,1}},QQ)
Hinv=inverse H

D=Hinv*M*H
H*M*Hinv

Hinv**Hinv

flattQ=time (Hinv**Hinv)*flattq*transpose(Hinv**Hinv);
flattQ



---------------------------------
restart
R=QQ[p_1,p_2,p_3,p_4,t]
--I=trim ideal{(p_1+p_2+p_3+p_4)*t,1-(p_1+p_2+p_3+p_4)}

--M=matrix{{1-(p_2+p_3+p_4)*t,p_2*t,p_3*t,p_4*t},{p_1*t,1-(p_1+p_3+p_4)*t,p_3*t,p_4*t},{p_1*t,p_2*t,1-(p_1+p_2+p_4)*t,p_4*t},{p_1*t,p_2*t,p_3*t,1-(p_1+p_2+p_3)*t}}
M=matrix{{1-(1-p_1)*t,p_2*t,p_3*t,p_4*t},{p_1*t,1-(1-p_2)*t,p_3*t,p_4*t},{p_1*t,p_2*t,1-(1-p_3)*t,p_4*t},{p_1*t,p_2*t,p_3*t,1-(1-p_4)*t}}

flattq=mutableMatrix(R,16,16)

for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=p_(i%4+1)*M_(i%4,j%4));      
    )

flattq=matrix(flattq)

--Base Angelica
H=transpose(matrix{{1,1,1,1},{-p_2,p_1,0,0},{-p_3,0,p_1,0},{-p_4,0,0,p_1}});

--K=frac(QQ[p_1,p_2,p_3,p_4,t]);
--H=sub(H,K);
--M=sub(M,K);
--Hinv=inverse(H);
--D=Hinv*M*H
--H*M*Hinv

flattQ=time (transpose(H)**transpose(H))*flattq*(H**H);
flattQ
netList (flatten entries flattQ)

netList toList apply(0..15,i->toList apply(0..15,j->(i,j,flattQ_(i,j))))
for i to 15 do (for j to 15 do (if flattQ_(i,j)==0 then print (i,j)))

extractFactors=factors->(
for i from 0 to #factors-1 list (toList factors#i)_0
);

for i to 15 do (for j to 15 do (if flattQ_(i,j)!=0 then print (i,j,extractFactors factor flattQ_(i,j))))

for j to 15 do (if flattQ_(0,j)!=0 then print (0,j,extractFactors factor flattQ_(0,j)))

for i to 15 do (for j to 15 do print (i,j,degree flattQ_(i,j)))

factor flattQ_(1,4)
factor flattQ_(12,2)
--
RQ=R/(1-(p_1+p_2+p_3+p_4))
flattQRQ=sub(flattQ,RQ);
for j to 15 do print (0,j,flattQRQ_(0,j))
for j to 15 do print (0,j,flattQ_(0,j))

sub(flattQ,{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})
sub(M,{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})

    
--Base Marta
restart
K=frac(QQ[p_1,p_2,p_3,p_4,t]);
M=matrix{{1-(1-p_1)*t,p_2*t,p_3*t,p_4*t},{p_1*t,1-(1-p_2)*t,p_3*t,p_4*t},{p_1*t,p_2*t,1-(1-p_3)*t,p_4*t},{p_1*t,p_2*t,p_3*t,1-(1-p_4)*t}}
--H=transpose(matrix{{1,1,1,1},{1/p_1,1/p_2,-1/p_3,-1/p_4},{1/p_1,-1/p_2,1/p_3,-1/p_4},{1/p_1,-1/p_2,-1/p_3,1/p_4}})
H=transpose(matrix{{4,4,4,4},{1/p_1,1/p_2,-1/p_3,-1/p_4},{1/p_1,-1/p_2,1/p_3,-1/p_4},{1/p_1,-1/p_2,-1/p_3,1/p_4}})

Hinv=inverse H
D=Hinv*M*H
H*M*Hinv


flattq=mutableMatrix(K,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=p_(i%4+1)*M_(i%4,j%4));      
    )
flattq=matrix(flattq)
flattQ=time (transpose(H)**transpose(H))*flattq*(H**H);


aux=flattQ
auxcol=submatrix(aux,{0,5,10,15})|submatrix(aux,{1,4,11,14})|submatrix(aux,{2,7,8,13})|submatrix(aux,{3,6,9,12})
auxrow=auxcol^{0,5,10,15}||auxcol^{1,4,11,14}||auxcol^{2,7,8,13}||auxcol^{3,6,9,12}
sub(auxrow,{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})

B11=submatrix(auxrow,{0,1,2,3},{0,1,2,3})
rank oo

B1=auxrow^{0,1,2,3}
rank B1 

B1^{0}==B1^{1}
B1^{1}==B1^{2}
B1^{1}==B1^{3}
--1,2,3 equal
auxrow^{4}==auxrow^{5} --true
auxrow^{4}==auxrow^{6} --false
auxrow^{4}==auxrow^{7} --false
auxrow^{6}==auxrow^{7} --true
--4,5 equal, 6,7 equal
auxrow^{8}==auxrow^{9} --false
auxrow^{8}==auxrow^{10} --true
auxrow^{9}==auxrow^{11} --true

auxrow^{12}==auxrow^{13} --false
auxrow^{12}==auxrow^{14} --false
auxrow^{12}==auxrow^{15} --true
auxrow^{13}==auxrow^{14} --true

auxrow_{0}==auxrow_{1} --false
auxrow_{1}==auxrow_{2} --true
auxrow_{1}==auxrow_{3} --true
--1,2,3 equal
auxrow_{4}==auxrow_{5} --true
auxrow_{4}==auxrow_{6} --false
auxrow_{4}==auxrow_{7} --false
auxrow_{6}==auxrow_{7} --true
--4,5 equal, 6,7 equal
auxrow_{8}==auxrow_{9} --false
auxrow_{8}==auxrow_{10} --true
auxrow_{9}==auxrow_{11} --true

auxrow_{12}==auxrow_{13} --false
auxrow_{12}==auxrow_{14} --false
auxrow_{12}==auxrow_{15} --true
auxrow_{13}==auxrow_{14} --true

gens K
R2=frac(QQ[gens K,a_1,a_2,a_3,a_4])
M1=matrix{{1,0,0,0},{0,1-a_1,0,0},{0,0,1-a_1,0},{0,0,0,1-a_1}}
M2=matrix{{1,0,0,0},{0,1-a_2,0,0},{0,0,1-a_2,0},{0,0,0,1-a_2}}
M3=matrix{{1,0,0,0},{0,1-a_3,0,0},{0,0,1-a_3,0},{0,0,0,1-a_3}}
M4=matrix{{1,0,0,0},{0,1-a_4,0,0},{0,0,1-a_4,0},{0,0,0,1-a_4}}
flattQ=sub(flattQ,R2);
flattP=transpose(M1**M2)*flattQ*(M3**M4);
auxP=flattP;
auxcolP=submatrix(auxP,{0,5,10,15})|submatrix(auxP,{1,4,11,14})|submatrix(auxP,{2,7,8,13})|submatrix(auxP,{3,6,9,12});
auxrowP=auxcolP^{0,5,10,15}||auxcolP^{1,4,11,14}||auxcolP^{2,7,8,13}||auxcolP^{3,6,9,12};

rank auxrowP^{0,1,2,3}
rank auxrowP_{0,1,2,3}

B12=submatrix(auxrow,{4,5,6,7},{0,1,2,3})
rank oo

sub(flattQ,{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})
sub(M,{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})


for j to 15 do (if flattQ_(0,j)!=0 then print (0,j,factor flattQ_(0,j)))

netList toList apply(flatten entries flattQ,i->factor numerator i)
netList toList apply(flatten entries flattQ,i->denominator i)
for i to 15 do (for j to 15 do (if flattQ_(i,j)==0 then print (i,j)))

apply(flatten entries flattQ,i->denominator i)
flattQR=sub(flattQ,R);


L=flatten entries flattQR;
Lsim=apply(flatten entries flattQR,i->(trim ideal{i,1-(p_1+p_2+p_3+p_4)})_1);

L_17
Lsim_17
aux=trim ideal{L_17,1-(p_1+p_2+p_3+p_4)}
netList aux_*

-----------
restart
K=frac(QQ[p_1,p_2,p_3,p_4,t]);
M=matrix{{1-(1-p_1)*t,p_2*t,p_3*t,p_4*t},{p_1*t,1-(1-p_2)*t,p_3*t,p_4*t},{p_1*t,p_2*t,1-(1-p_3)*t,p_4*t},{p_1*t,p_2*t,p_3*t,1-(1-p_4)*t}}
--H=transpose(matrix{{1,1,1,1},{1/p_1,1/p_2,-1/p_3,-1/p_4},{1/p_1,-1/p_2,1/p_3,-1/p_4},{1/p_1,-1/p_2,-1/p_3,1/p_4}})
H=transpose(matrix{{1/p_1,1/p_1,1/p_1,1/p_1},{1/p_1,1/p_2,-1/p_3,-1/p_4},{1/p_1,-1/p_2,1/p_3,-1/p_4},{1/p_1,-1/p_2,-1/p_3,1/p_4}})

Hinv=inverse H
D=Hinv*M*H
H*M*Hinv

aux=flattQ;


auxcol=submatrix(aux,{0,5,10,15})|submatrix(aux,{1,4,11,14})|submatrix(aux,{2,7,8,13})|submatrix(aux,{3,6,9,12})
auxrow=auxcol^{0,5,10,15}||auxcol^{1,4,11,14}||auxcol^{2,7,8,13}||auxcol^{3,6,9,12}
sub(auxrow,{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})

submatrix(auxrow,{0,1,2,3},{0,1,2,3})

flattq=mutableMatrix(K,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=p_(i%4+1)*M_(i%4,j%4));      
    )
flattq=matrix(flattq)
flattQ=time (transpose(H)**transpose(H))*flattq*(H**H);

sub(flattQ,{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})
sub(M,{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})


for j to 15 do (if flattQ_(0,j)!=0 then print (0,j,factor flattQ_(0,j)))
for i to 15 do (for j to 15 do (if flattQ_(i,j)==0 then print (i,j)))
