% Plotting the ``bridge'' truss structure
bridge; % Initialize data fields
figure('name','bridge');
trussplot(pos,top);
print -depsc2 '../PICTURES/bridgetruss.eps';
