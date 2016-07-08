function Uout = circsim(U_in,alpha,beta)
    Uout = [];
    R    = 1;
    U_t  = 0.5;
    
    for m=1:length(U_in)
        U=rand(3,1);
        h=[1 1 1]';
    
        while max(abs(h)) > 1.0e-6*norm(U)
            J = 1/R*[  3 -1 -1;
                      -1  3 -1;
     -1 -1  3+(R*alpha*beta)/U_t*exp(beta*(U(3)-U_in(m))/U_t) ];
              
            f = [ (3*U(1)-U(2)-U(3))/R;
                  (3*U(2)-U(1)-U(3)-U_in(m))/R;
  (3*U(3)-U(1)-U(2))/R + alpha*(exp(beta*(U(3)-U_in(m))/U_t)-1) ];
        
            h = J\f;
            U = U-h;
        end
        Uout = [Uout U(1)];
    end