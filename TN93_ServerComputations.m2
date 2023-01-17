K=frac(QQ[p_A,p_C,p_G,p_T]);
gensp={p_((A,A),(A,A)), p_((A,A),(A,T)), p_((A,A),(T,A)), p_((A,A),(T,G)), p_((A,A),(G,T)), p_((A,A),(T,T)), p_((A,A),(G,G)), p_((A,A),(A,G)), p_((A,A),(G,A)), p_((A,A),(C,C)),
      p_((A,A),(A,C)), p_((A,A),(C,A)), p_((A,A),(C,G)), p_((A,A),(G,C)), p_((A,A),(T,C)), p_((A,A),(C,T)), p_((A,T),(A,A)), p_((A,T),(A,T)), p_((A,T),(T,A)), p_((A,T),(T,G)),
      p_((A,T),(G,T)), p_((A,T),(T,T)), p_((A,T),(G,G)), p_((A,T),(A,G)), p_((A,T),(G,A)), p_((A,T),(C,C)), p_((A,T),(A,C)), p_((A,T),(C,A)), p_((A,T),(C,G)), p_((A,T),(G,C)),
      p_((A,T),(T,C)), p_((A,T),(C,T)), p_((T,A),(A,A)), p_((T,A),(A,T)), p_((T,A),(T,A)), p_((T,A),(T,G)), p_((T,A),(G,T)), p_((T,A),(T,T)), p_((T,A),(G,G)), p_((T,A),(A,G)),
      p_((T,A),(G,A)), p_((T,A),(C,C)), p_((T,A),(A,C)), p_((T,A),(C,A)), p_((T,A),(C,G)), p_((T,A),(G,C)), p_((T,A),(T,C)), p_((T,A),(C,T)), p_((T,G),(A,A)), p_((T,G),(A,T)),
      p_((T,G),(T,A)), p_((T,G),(T,G)), p_((T,G),(G,T)), p_((T,G),(T,T)), p_((T,G),(G,G)), p_((T,G),(A,G)), p_((T,G),(G,A)), p_((T,G),(C,C)), p_((T,G),(A,C)), p_((T,G),(C,A)),
      p_((T,G),(C,G)), p_((T,G),(G,C)), p_((T,G),(T,C)), p_((T,G),(C,T)), p_((G,T),(A,A)), p_((G,T),(A,T)), p_((G,T),(T,A)), p_((G,T),(T,G)), p_((G,T),(G,T)), p_((G,T),(T,T)),
      p_((G,T),(G,G)), p_((G,T),(A,G)), p_((G,T),(G,A)), p_((G,T),(C,C)), p_((G,T),(A,C)), p_((G,T),(C,A)), p_((G,T),(C,G)), p_((G,T),(G,C)), p_((G,T),(T,C)), p_((G,T),(C,T)),
      p_((T,T),(A,A)), p_((T,T),(A,T)), p_((T,T),(T,A)), p_((T,T),(T,G)), p_((T,T),(G,T)), p_((T,T),(T,T)), p_((T,T),(G,G)), p_((T,T),(A,G)), p_((T,T),(G,A)), p_((T,T),(C,C)),
      p_((T,T),(A,C)), p_((T,T),(C,A)), p_((T,T),(C,G)), p_((T,T),(G,C)), p_((T,T),(T,C)), p_((T,T),(C,T)), p_((G,G),(A,A)), p_((G,G),(A,T)), p_((G,G),(T,A)), p_((G,G),(T,G)),
      p_((G,G),(G,T)), p_((G,G),(T,T)), p_((G,G),(G,G)), p_((G,G),(A,G)), p_((G,G),(G,A)), p_((G,G),(C,C)), p_((G,G),(A,C)), p_((G,G),(C,A)), p_((G,G),(C,G)), p_((G,G),(G,C)),
      p_((G,G),(T,C)), p_((G,G),(C,T)), p_((A,G),(A,A)), p_((A,G),(A,T)), p_((A,G),(T,A)), p_((A,G),(T,G)), p_((A,G),(G,T)), p_((A,G),(T,T)), p_((A,G),(G,G)), p_((A,G),(A,G)),
      p_((A,G),(G,A)), p_((A,G),(C,C)), p_((A,G),(A,C)), p_((A,G),(C,A)), p_((A,G),(C,G)), p_((A,G),(G,C)), p_((A,G),(T,C)), p_((A,G),(C,T)), p_((G,A),(A,A)), p_((G,A),(A,T)),
      p_((G,A),(T,A)), p_((G,A),(T,G)), p_((G,A),(G,T)), p_((G,A),(T,T)), p_((G,A),(G,G)), p_((G,A),(A,G)), p_((G,A),(G,A)), p_((G,A),(C,C)), p_((G,A),(A,C)), p_((G,A),(C,A)),
      p_((G,A),(C,G)), p_((G,A),(G,C)), p_((G,A),(T,C)), p_((G,A),(C,T)), p_((C,C),(A,A)), p_((C,C),(A,T)), p_((C,C),(T,A)), p_((C,C),(T,G)), p_((C,C),(G,T)), p_((C,C),(T,T)),
      p_((C,C),(G,G)), p_((C,C),(A,G)), p_((C,C),(G,A)), p_((C,C),(C,C)), p_((C,C),(A,C)), p_((C,C),(C,A)), p_((C,C),(C,G)), p_((C,C),(G,C)), p_((C,C),(T,C)), p_((C,C),(C,T)),
      p_((A,C),(A,A)), p_((A,C),(A,T)), p_((A,C),(T,A)), p_((A,C),(T,G)), p_((A,C),(G,T)), p_((A,C),(T,T)), p_((A,C),(G,G)), p_((A,C),(A,G)), p_((A,C),(G,A)), p_((A,C),(C,C)),
      p_((A,C),(A,C)), p_((A,C),(C,A)), p_((A,C),(C,G)), p_((A,C),(G,C)), p_((A,C),(T,C)), p_((A,C),(C,T)), p_((C,A),(A,A)), p_((C,A),(A,T)), p_((C,A),(T,A)), p_((C,A),(T,G)),
      p_((C,A),(G,T)), p_((C,A),(T,T)), p_((C,A),(G,G)), p_((C,A),(A,G)), p_((C,A),(G,A)), p_((C,A),(C,C)), p_((C,A),(A,C)), p_((C,A),(C,A)), p_((C,A),(C,G)), p_((C,A),(G,C)),
      p_((C,A),(T,C)), p_((C,A),(C,T)), p_((C,G),(A,A)), p_((C,G),(A,T)), p_((C,G),(T,A)), p_((C,G),(T,G)), p_((C,G),(G,T)), p_((C,G),(T,T)), p_((C,G),(G,G)), p_((C,G),(A,G)),
      p_((C,G),(G,A)), p_((C,G),(C,C)), p_((C,G),(A,C)), p_((C,G),(C,A)), p_((C,G),(C,G)), p_((C,G),(G,C)), p_((C,G),(T,C)), p_((C,G),(C,T)), p_((G,C),(A,A)), p_((G,C),(A,T)),
      p_((G,C),(T,A)), p_((G,C),(T,G)), p_((G,C),(G,T)), p_((G,C),(T,T)), p_((G,C),(G,G)), p_((G,C),(A,G)), p_((G,C),(G,A)), p_((G,C),(C,C)), p_((G,C),(A,C)), p_((G,C),(C,A)),
      p_((G,C),(C,G)), p_((G,C),(G,C)), p_((G,C),(T,C)), p_((G,C),(C,T)), p_((T,C),(A,A)), p_((T,C),(A,T)), p_((T,C),(T,A)), p_((T,C),(T,G)), p_((T,C),(G,T)), p_((T,C),(T,T)),
      p_((T,C),(G,G)), p_((T,C),(A,G)), p_((T,C),(G,A)), p_((T,C),(C,C)), p_((T,C),(A,C)), p_((T,C),(C,A)), p_((T,C),(C,G)), p_((T,C),(G,C)), p_((T,C),(T,C)), p_((T,C),(C,T)),
      p_((C,T),(A,A)), p_((C,T),(A,T)), p_((C,T),(T,A)), p_((C,T),(T,G)), p_((C,T),(G,T)), p_((C,T),(T,T)), p_((C,T),(G,G)), p_((C,T),(A,G)), p_((C,T),(G,A)), p_((C,T),(C,C)),
      p_((C,T),(A,C)), p_((C,T),(C,A)), p_((C,T),(C,G)), p_((C,T),(G,C)), p_((C,T),(T,C)), p_((C,T),(C,T))};
gensl={l_(1,A), l_(1,C), l_(1,G), l_(1,T), l_(2,A), l_(2,C), l_(2,G), l_(2,T), l_(3,A),
      l_(3,C), l_(3,G), l_(3,T), l_(4,A), l_(4,C), l_(4,G), l_(4,T), l_(5,A), l_(5,C), l_(5,G), l_(5,T)};
deg={5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
R=K[gensp,gensl,Degrees=>deg]

I=ideal(p_((C,T),(C,T)),p_((C,T),(T,C)),p_((C,T),(G,C)),p_((C,T),(C,G)),p_((C,T),(C,A)),p_((C,T),(A,C)),p_((C,T),(C,C)),p_((C,T),(G,A)),p_((C,T),(A,G)),p_((C,T),(G,G)),p_((C,T),(T,T)),
p_((C,T),(G,T)),p_((C,T),(T,G)),p_((C,T),(T,A)),p_((C,T),(A,T)),p_((C,T),(A,A)),p_((T,C),(C,T)),p_((T,C),(T,C)),p_((T,C),(G,C)),p_((T,C),(C,G)),p_((T,C),(C,A)),p_((T,C),(A,C)),p_((T,C),(C,C)),
p_((T,C),(G,A)),p_((T,C),(A,G)),p_((T,C),(G,G)),p_((T,C),(T,T)),p_((T,C),(G,T)),p_((T,C),(T,G)),p_((T,C),(T,A)),p_((T,C),(A,T)),p_((T,C),(A,A)),p_((G,C),(C,T)),p_((G,C),(T,C)),
p_((G,C),(G,C))+((-1)/(p_C^2*p_T+p_C*p_T^2))*l_(1,G)*l_(2,C)*l_(3,G)*l_(4,C)*l_(5,C),p_((G,C),(C,G))+((-1)/(p_C^2*p_T+p_C*p_T^2))*l_(1,G)*l_(2,C)*l_(3,C)*l_(4,G)*l_(5,C),
p_((G,C),(C,A))+(4/(p_C*p_T))*l_(1,G)*l_(2,C)*l_(3,C)*l_(4,A)*l_(5,C),p_((G,C),(A,C))+(4/(p_C*p_T))*l_(1,G)*l_(2,C)*l_(3,A)*l_(4,C)*l_(5,C),
p_((G,C),(C,C))+((-p_C+p_T)/(p_C^2*p_T^2))*l_(1,G)*l_(2,C)*l_(3,C)*l_(4,C)*l_(5,C),p_((G,C),(G,A)),p_((G,C),(A,G)),p_((G,C),(G,G)),p_((G,C),(T,T)),p_((G,C),(G,T)),p_((G,C),(T,G)),p_((G,C),(T,A)),p_((G,C),(A,T)),p_((G,C),(A,A)),p_((C,G),(C,T)),
p_((C,G),(T,C)),p_((C,G),(G,C))+((-1)/(p_C^2*p_T+p_C*p_T^2))*l_(1,C)*l_(2,G)*l_(3,G)*l_(4,C)*l_(5,C),p_((C,G),(C,G))+((-1)/(p_C^2*p_T+p_C*p_T^2))*l_(1,C)*l_(2,G)*l_(3,C)*l_(4,G)*l_(5,C),
p_((C,G),(C,A))+(4/(p_C*p_T))*l_(1,C)*l_(2,G)*l_(3,C)*l_(4,A)*l_(5,C),p_((C,G),(A,C))+(4/(p_C*p_T))*l_(1,C)*l_(2,G)*l_(3,A)*l_(4,C)*l_(5,C),p_((C,G),(C,C))+((-p_C+p_T)/(p_C^2*p_T^2))*l_(1,C)*l_(2,G)*l_(3,C)*l_(4,C)*l_(5,C),p_((C,G),(G,A)),p_((C,G),(A,G)),p_((C,G),(G,G)),p_((C,G),(T,T)),p_((C,G),(G,T)),p_((C,G),(T,G)),p_((C,G),(T,A)),p_((C,G),(A,T)),
p_((C,G),(A,A)),p_((C,A),(C,T)),p_((C,A),(T,C)),p_((C,A),(G,C))+(4/(p_C*p_T))*l_(1,C)*l_(2,A)*l_(3,G)*l_(4,C)*l_(5,C),p_((C,A),(C,G))+(4/(p_C*p_T))*l_(1,C)*l_(2,A)*l_(3,C)*l_(4,G)*l_(5,C),
p_((C,A),(C,A))+((-16*p_C-16*p_T)/(p_C*p_T))*l_(1,C)*l_(2,A)*l_(3,C)*l_(4,A)*l_(5,C),p_((C,A),(A,C))+((-16*p_C-16*p_T)/(p_C*p_T))*l_(1,C)*l_(2,A)*l_(3,A)*l_(4,C)*l_(5,C),
p_((C,A),(C,C))+((4*p_C^2-4*p_T^2)/(p_C^2*p_T^2))*l_(1,C)*l_(2,A)*l_(3,C)*l_(4,C)*l_(5,C),p_((C,A),(G,A)),p_((C,A),(A,G)),p_((C,A),(G,G)),p_((C,A),(T,T)),p_((C,A),(G,T)),p_((C,A),(T,G)),
p_((C,A),(T,A)),p_((C,A),(A,T)),p_((C,A),(A,A)),p_((A,C),(C,T)),p_((A,C),(T,C)),p_((A,C),(G,C))+(4/(p_C*p_T))*l_(1,A)*l_(2,C)*l_(3,G)*l_(4,C)*l_(5,C),
p_((A,C),(C,G))+(4/(p_C*p_T))*l_(1,A)*l_(2,C)*l_(3,C)*l_(4,G)*l_(5,C),p_((A,C),(C,A))+((-16*p_C-16*p_T)/(p_C*p_T))*l_(1,A)*l_(2,C)*l_(3,C)*l_(4,A)*l_(5,C),
p_((A,C),(A,C))+((-16*p_C-16*p_T)/(p_C*p_T))*l_(1,A)*l_(2,C)*l_(3,A)*l_(4,C)*l_(5,C),p_((A,C),(C,C))+((4*p_C^2-4*p_T^2)/(p_C^2*p_T^2))*l_(1,A)*l_(2,C)*l_(3,C)*l_(4,C)*l_(5,C),p_((A,C),(G,A)),p_((A,C),(A,G)),p_((A,C),(G,G)),p_((A,C),(T,T)),
p_((A,C),(G,T)),p_((A,C),(T,G)),p_((A,C),(T,A)),p_((A,C),(A,T)),p_((A,C),(A,A)),p_((C,C),(C,T)),p_((C,C),(T,C)),
p_((C,C),(G,C))+((-p_C+p_T)/(p_C^2*p_T^2))*l_(1,C)*l_(2,C)*l_(3,G)*l_(4,C)*l_(5,C),p_((C,C),(C,G))+((-p_C+p_T)/(p_C^2*p_T^2))*l_(1,C)*l_(2,C)*l_(3,C)*l_(4,G)*l_(5,C),p_((C,C),(C,A))+((4*p_C^2-4*p_T^2)/(p_C^2*p_T^2))*l_(1,C)*l_(2,C)*l_(3,C)*l_(4,A)*l_(5,C),
p_((C,C),(A,C))+((4*p_C^2-4*p_T^2)/(p_C^2*p_T^2))*l_(1,C)*l_(2,C)*l_(3,A)*l_(4,C)*l_(5,C),p_((C,C),(C,C))+((-p_C^2-2*p_C*p_T-p_T^2)/(p_A*p_C^2*p_T^2+p_C^3*p_T^2+p_C^2*p_G*p_T^2+p_C^2*p_T^3))*l_(1,C)*l_(2,C)*l_(3,C)*l_(4,C)*l_(5,A)+((-p_C^3+p_C^2*p_T+p_C*p_T^2-p_T^3)/(p_C^3*p_T^3))*l_(1,C)*l_(2,C)*l_(3,C)*l_(4,C)*l_(5,C)+((-p_A*p_C-p_C*p_G-p_A*p_T-p_G*p_T)/(p_A*p_C^2*p_T^2+p_C^3*p_T^2+p_C^2*p_G*p_T^2+p_C^2*p_T^3))*l_(1,C)*l_(2,C)*l_(3,C)*l_(4,C)*l_(5,G),
p_((C,C),(G,A))+(4/(p_C*p_T))*l_(1,C)*l_(2,C)*l_(3,G)*l_(4,A)*l_(5,G),p_((C,C),(A,G))+(4/(p_C*p_T))*l_(1,C)*l_(2,C)*l_(3,A)*l_(4,G)*l_(5,G),p_((C,C),(G,G))+((-1)/(p_A*p_C*p_T+p_C*p_G*p_T))*l_(1,C)*l_(2,C)*l_(3,G)*l_(4,G)*l_(5,A)+((-p_A+p_C-p_G+p_T)/(p_A*p_C^2*p_T+p_C^2*p_G*p_T+p_A*p_C*p_T^2+p_C*p_G*p_T^2))*l_(1,C)*l_(2,C)*l_(3,G)*l_(4,G)*l_(5,G),
p_((C,C),(T,T))+((-p_A*p_C-p_C*p_G-p_A*p_T-p_G*p_T)/(p_A^2*p_C*p_G*p_T+p_A*p_C^2*p_G*p_T+p_A*p_C*p_G^2*p_T+p_A*p_C*p_G*p_T^2))*l_(1,C)*l_(2,C)*l_(3,T)*l_(4,T)*l_(5,A)+((p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)/(p_A^2*p_C*p_G*p_T+p_A*p_C^2*p_G*p_T+p_A*p_C*p_G^2*p_T+p_A*p_C*p_G*p_T^2))*l_(1,C)*l_(2,C)*l_(3,T)*l_(4,T)*l_(5,G),
p_((C,C),(G,T)),p_((C,C),(T,G)),p_((C,C),(T,A)),p_((C,C),(A,T)),p_((C,C),(A,A))+((-16*p_C-16*p_T)/(p_C*p_T))*l_(1,C)*l_(2,C)*l_(3,A)*l_(4,A)*l_(5,A),p_((G,A),(C,T)),p_((G,A),(T,C)),p_((G,A),(G,C)),
p_((G,A),(C,G)),p_((G,A),(C,A)),p_((G,A),(A,C)),p_((G,A),(C,C))+(4/(p_C*p_T))*l_(1,G)*l_(2,A)*l_(3,C)*l_(4,C)*l_(5,G),p_((G,A),(G,A))+((-16*p_A-16*p_C-16*p_G-16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(1,G)*l_(2,A)*l_(3,G)*l_(4,A)*l_(5,G),p_((G,A),(A,G))+((-16*p_A-16*p_C-16*p_G-16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(1,G)*l_(2,A)*l_(3,A)*l_(4,G)*l_(5,G),
p_((G,A),(G,G))+((4*p_A^2-4*p_C^2+8*p_A*p_G+4*p_G^2-8*p_C*p_T-4*p_T^2)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(1,G)*l_(2,A)*l_(3,G)*l_(4,G)*l_(5,G),p_((G,A),(T,T))+((-4)/(p_A*p_G))*l_(1,G)*l_(2,A)*l_(3,T)*l_(4,T)*l_(5,G),
p_((G,A),(G,T)),p_((G,A),(T,G)),p_((G,A),(T,A)),p_((G,A),(A,T)),p_((G,A),(A,A)),p_((A,G),(C,T)),p_((A,G),(T,C)),p_((A,G),(G,C)),p_((A,G),(C,G)),p_((A,G),(C,A)),p_((A,G),(A,C)),
p_((A,G),(C,C))+(4/(p_C*p_T))*l_(1,A)*l_(2,G)*l_(3,C)*l_(4,C)*l_(5,G),p_((A,G),(G,A))+((-16*p_A-16*p_C-16*p_G-16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(1,A)*l_(2,G)*l_(3,G)*l_(4,A)*l_(5,G),
p_((A,G),(A,G))+((-16*p_A-16*p_C-16*p_G-16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(1,A)*l_(2,G)*l_(3,A)*l_(4,G)*l_(5,G),p_((A,G),(G,G))+((4*p_A^2-4*p_C^2+8*p_A*p_G+4*p_G^2-8*p_C*p_T-4*p_T^2)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(1,A)*l_(2,G)*l_(3,G)*l_(4,G)*l_(5,G),
p_((A,G),(T,T))+((-4)/(p_A*p_G))*l_(1,A)*l_(2,G)*l_(3,T)*l_(4,T)*l_(5,G),p_((A,G),(G,T)),p_((A,G),(T,G)),p_((A,G),(T,A)),p_((A,G),(A,T)),p_((A,G),(A,A)),p_((G,G),(C,T)),
p_((G,G),(T,C)),p_((G,G),(G,C)),p_((G,G),(C,G)),p_((G,G),(C,A)),p_((G,G),(A,C)),p_((G,G),(C,C))+((-1)/(p_A*p_C*p_T+p_C*p_G*p_T))*l_(1,G)*l_(2,G)*l_(3,C)*l_(4,C)*l_(5,A)+((-p_A+p_C-p_G+p_T)/(p_A*p_C^2*p_T+p_C^2*p_G*p_T+p_A*p_C*p_T^2+p_C*p_G*p_T^2))*l_(1,G)*l_(2,G)*l_(3,C)*l_(4,C)*l_(5,G),
p_((G,G),(G,A))+((4*p_A^2-4*p_C^2+8*p_A*p_G+4*p_G^2-8*p_C*p_T-4*p_T^2)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(1,G)*l_(2,G)*l_(3,G)*l_(4,A)*l_(5,G),p_((G,G),(A,G))+((4*p_A^2-4*p_C^2+8*p_A*p_G+4*p_G^2-8*p_C*p_T-4*p_T^2)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(1,G)*l_(2,G)*l_(3,A)*l_(4,G)*l_(5,G),
p_((G,G),(G,G))+((-p_A-p_C-p_G-p_T)/(p_A^2*p_C^2+2*p_A*p_C^2*p_G+p_C^2*p_G^2+2*p_A^2*p_C*p_T+4*p_A*p_C*p_G*p_T+2*p_C*p_G^2*p_T+p_A^2*p_T^2+2*p_A*p_G*p_T^2+p_G^2*p_T^2))*l_(1,G)*l_(2,G)*l_(3,G)*l_(4,G)*l_(5,A)+((-p_A^3+p_A^2*p_C+p_A*p_C^2-p_C^3-3*p_A^2*p_G+2*p_A*p_C*p_G+p_C^2*p_G-3*p_A*p_G^2+p_C*p_G^2-p_G^3+p_A^2*p_T+2*p_A*p_C*p_T-3*p_C^2*p_T+2*p_A*p_G*p_T+2*p_C*p_G*p_T+p_G^2*p_T+p_A*p_T^2-3*p_C*p_T^2+p_G*p_T^2-p_T^3)/(p_A^3*p_C^3+3*p_A^2*p_C^3*p_G+3*p_A*p_C^3*p_G^2+p_C^3*p_G^3+3*p_A^3*p_C^2*p_T+9*p_A^2*p_C^2*p_G*p_T+9*p_A*p_C^2*p_G^2*p_T+3*p_C^2*p_G^3*p_T+3*p_A^3*p_C*p_T^2+9*p_A^2*p_C*p_G*p_T^2+9*p_A*p_C*p_G^2*p_T^2+3*p_C*p_G^3*p_T^2+p_A^3*p_T^3+3*p_A^2*p_G*p_T^3+3*p_A*p_G^2*p_T^3+p_G^3*p_T^3))*l_(1,G)*l_(2,G)*l_(3,G)*l_(4,G)*l_(5,G),
p_((G,G),(T,T))+((-1)/(p_A*p_C*p_G+p_A*p_G*p_T))*l_(1,G)*l_(2,G)*l_(3,T)*l_(4,T)*l_(5,A)+((p_A-p_C+p_G-p_T)/(p_A^2*p_C*p_G+p_A*p_C*p_G^2+p_A^2*p_G*p_T+p_A*p_G^2*p_T))*l_(1,G)*l_(2,G)*l_(3,T)*l_(4,T)*l_(5,G),p_((G,G),(G,T)),
p_((G,G),(T,G)),p_((G,G),(T,A)),p_((G,G),(A,T)),p_((G,G),(A,A))+((-16*p_A-16*p_C-16*p_G-16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(1,G)*l_(2,G)*l_(3,A)*l_(4,A)*l_(5,A),
p_((T,T),(C,T)),p_((T,T),(T,C)),p_((T,T),(G,C)),p_((T,T),(C,G)),p_((T,T),(C,A)),p_((T,T),(A,C)),p_((T,T),(C,C))+((-p_A*p_C-p_C*p_G-p_A*p_T-p_G*p_T)/(p_A^2*p_C*p_G*p_T+p_A*p_C^2*p_G*p_T+p_A*p_C*p_G^2*p_T+p_A*p_C*p_G*p_T^2))*l_(1,T)*l_(2,T)*l_(3,C)*l_(4,C)*l_(5,A)+((p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T)/(p_A^2*p_C*p_G*p_T+p_A*p_C^2*p_G*p_T+p_A*p_C*p_G^2*p_T+p_A*p_C*p_G*p_T^2))*l_(1,T)*l_(2,T)*l_(3,C)*l_(4,C)*l_(5,G),
p_((T,T),(G,A))+((-4)/(p_A*p_G))*l_(1,T)*l_(2,T)*l_(3,G)*l_(4,A)*l_(5,G),p_((T,T),(A,G))+((-4)/(p_A*p_G))*l_(1,T)*l_(2,T)*l_(3,A)*l_(4,G)*l_(5,G),
p_((T,T),(G,G))+((-1)/(p_A*p_C*p_G+p_A*p_G*p_T))*l_(1,T)*l_(2,T)*l_(3,G)*l_(4,G)*l_(5,A)+((p_A-p_C+p_G-p_T)/(p_A^2*p_C*p_G+p_A*p_C*p_G^2+p_A^2*p_G*p_T+p_A*p_G^2*p_T))*l_(1,T)*l_(2,T)*l_(3,G)*l_(4,G)*l_(5,G),
p_((T,T),(T,T))+((-p_A^2-2*p_A*p_G-p_G^2)/(p_A^3*p_G^2+p_A^2*p_C*p_G^2+p_A^2*p_G^3+p_A^2*p_G^2*p_T))*l_(1,T)*l_(2,T)*l_(3,T)*l_(4,T)*l_(5,A)+((-p_A*p_C-p_C*p_G-p_A*p_T-p_G*p_T)/(p_A^3*p_G^2+p_A^2*p_C*p_G^2+p_A^2*p_G^3+p_A^2*p_G^2*p_T))*l_(1,T)*l_(2,T)*l_(3,T)*l_(4,T)*l_(5,G)+((-p_A^3+p_A^2*p_G+p_A*p_G^2-p_G^3)/(p_A^3*p_G^3))*l_(1,T)*l_(2,T)*l_(3,T)*l_(4,T)*l_(5,T),
p_((T,T),(G,T))+((p_A-p_G)/(p_A^2*p_G^2))*l_(1,T)*l_(2,T)*l_(3,G)*l_(4,T)*l_(5,T),p_((T,T),(T,G))+((p_A-p_G)/(p_A^2*p_G^2))*l_(1,T)*l_(2,T)*l_(3,T)*l_(4,G)*l_(5,T),
p_((T,T),(T,A))+((4*p_A^2-4*p_G^2)/(p_A^2*p_G^2))*l_(1,T)*l_(2,T)*l_(3,T)*l_(4,A)*l_(5,T),p_((T,T),(A,T))+((4*p_A^2-4*p_G^2)/(p_A^2*p_G^2))*l_(1,T)*l_(2,T)*l_(3,A)*l_(4,T)*l_(5,T),
p_((T,T),(A,A))+((-16*p_A-16*p_G)/(p_A*p_G))*l_(1,T)*l_(2,T)*l_(3,A)*l_(4,A)*l_(5,A),p_((G,T),(C,T)),p_((G,T),(T,C)),p_((G,T),(G,C)),p_((G,T),(C,G)),p_((G,T),(C,A)),p_((G,T),(A,C)),
p_((G,T),(C,C)),p_((G,T),(G,A)),p_((G,T),(A,G)),p_((G,T),(G,G)),p_((G,T),(T,T))+((p_A-p_G)/(p_A^2*p_G^2))*l_(1,G)*l_(2,T)*l_(3,T)*l_(4,T)*l_(5,T),
p_((G,T),(G,T))+((-1)/(p_A^2*p_G+p_A*p_G^2))*l_(1,G)*l_(2,T)*l_(3,G)*l_(4,T)*l_(5,T),p_((G,T),(T,G))+((-1)/(p_A^2*p_G+p_A*p_G^2))*l_(1,G)*l_(2,T)*l_(3,T)*l_(4,G)*l_(5,T),
p_((G,T),(T,A))+((-4)/(p_A*p_G))*l_(1,G)*l_(2,T)*l_(3,T)*l_(4,A)*l_(5,T),p_((G,T),(A,T))+((-4)/(p_A*p_G))*l_(1,G)*l_(2,T)*l_(3,A)*l_(4,T)*l_(5,T),p_((G,T),(A,A)),p_((T,G),(C,T)),
p_((T,G),(T,C)),p_((T,G),(G,C)),p_((T,G),(C,G)),p_((T,G),(C,A)),p_((T,G),(A,C)),p_((T,G),(C,C)),p_((T,G),(G,A)),p_((T,G),(A,G)),p_((T,G),(G,G)),p_((T,G),(T,T))+((p_A-p_G)/(p_A^2*p_G^2))*l_(1,T)*l_(2,G)*l_(3,T)*l_(4,T)*l_(5,T),p_((T,G),(G,T))+((-1)/(p_A^2*p_G+p_A*p_G^2))*l_(1,T)*l_(2,G)*l_(3,G)*l_(4,T)*l_(5,T),p_((T,G),(T,G))+((-1)/(p_A^2*p_G+p_A*p_G^2))*l_(1,T)*l_(2,G)*l_(3,T)*l_(4,G)*l_(5,T),
p_((T,G),(T,A))+((-4)/(p_A*p_G))*l_(1,T)*l_(2,G)*l_(3,T)*l_(4,A)*l_(5,T),p_((T,G),(A,T))+((-4)/(p_A*p_G))*l_(1,T)*l_(2,G)*l_(3,A)*l_(4,T)*l_(5,T),p_((T,G),(A,A)),p_((T,A),(C,T)),p_((T,A),(T,C)),
p_((T,A),(G,C)),p_((T,A),(C,G)),p_((T,A),(C,A)),p_((T,A),(A,C)),p_((T,A),(C,C)),p_((T,A),(G,A)),p_((T,A),(A,G)),p_((T,A),(G,G)),p_((T,A),(T,T))+((4*p_A^2-4*p_G^2)/(p_A^2*p_G^2))*l_(1,T)*l_(2,A)*l_(3,T)*l_(4,T)*l_(5,T),
p_((T,A),(G,T))+((-4)/(p_A*p_G))*l_(1,T)*l_(2,A)*l_(3,G)*l_(4,T)*l_(5,T),p_((T,A),(T,G))+((-4)/(p_A*p_G))*l_(1,T)*l_(2,A)*l_(3,T)*l_(4,G)*l_(5,T),p_((T,A),(T,A))+((-16*p_A-16*p_G)/(p_A*p_G))*l_(1,T)*l_(2,A)*l_(3,T)*l_(4,A)*l_(5,T),
p_((T,A),(A,T))+((-16*p_A-16*p_G)/(p_A*p_G))*l_(1,T)*l_(2,A)*l_(3,A)*l_(4,T)*l_(5,T),p_((T,A),(A,A)),p_((A,T),(C,T)),p_((A,T),(T,C)),p_((A,T),(G,C)),p_((A,T),(C,G)),p_((A,T),(C,A)),p_((A,T),(A,C)),p_((A,T),(C,C)),
p_((A,T),(G,A)),p_((A,T),(A,G)),p_((A,T),(G,G)),p_((A,T),(T,T))+((4*p_A^2-4*p_G^2)/(p_A^2*p_G^2))*l_(1,A)*l_(2,T)*l_(3,T)*l_(4,T)*l_(5,T),p_((A,T),(G,T))+((-4)/(p_A*p_G))*l_(1,A)*l_(2,T)*l_(3,G)*l_(4,T)*l_(5,T),
p_((A,T),(T,G))+((-4)/(p_A*p_G))*l_(1,A)*l_(2,T)*l_(3,T)*l_(4,G)*l_(5,T),p_((A,T),(T,A))+((-16*p_A-16*p_G)/(p_A*p_G))*l_(1,A)*l_(2,T)*l_(3,T)*l_(4,A)*l_(5,T),p_((A,T),(A,T))+((-16*p_A-16*p_G)/(p_A*p_G))*l_(1,A)*l_(2,T)*l_(3,A)*l_(4,T)*l_(5,T),
p_((A,T),(A,A)),p_((A,A),(C,T)),p_((A,A),(T,C)),p_((A,A),(G,C)),p_((A,A),(C,G)),p_((A,A),(C,A)),p_((A,A),(A,C)),p_((A,A),(C,C))+((-16*p_C-16*p_T)/(p_C*p_T))*l_(1,A)*l_(2,A)*l_(3,C)*l_(4,C)*l_(5,A),
p_((A,A),(G,A)),p_((A,A),(A,G)),p_((A,A),(G,G))+((-16*p_A-16*p_C-16*p_G-16*p_T)/(p_A*p_C+p_C*p_G+p_A*p_T+p_G*p_T))*l_(1,A)*l_(2,A)*l_(3,G)*l_(4,G)*l_(5,A),
p_((A,A),(T,T))+((-16*p_A-16*p_G)/(p_A*p_G))*l_(1,A)*l_(2,A)*l_(3,T)*l_(4,T)*l_(5,A),p_((A,A),(G,T)),p_((A,A),(T,G)),p_((A,A),(T,A)),p_((A,A),(A,T)),p_((A,A),(A,A))+(-256*p_A-256*p_C-256*p_G-256*p_T)*l_(1,A)*l_(2,A)*l_(3,A)*l_(4,A)*l_(5,A));

J=time eliminate(drop(gens R,256),I);

"EliminationIdeal.txt" << toString J << endl << close




