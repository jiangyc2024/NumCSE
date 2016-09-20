v = v/norm(v); 
w = w - dot(v,w)*v; 
w = w/norm(w);
