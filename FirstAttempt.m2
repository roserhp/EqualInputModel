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
