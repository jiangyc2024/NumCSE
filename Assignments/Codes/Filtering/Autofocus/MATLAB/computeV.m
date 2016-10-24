function V = computeV (Bh)
[n,m] = size(Bh);
[K1 K2] = meshgrid( circshift((1:m)-m/2,round(m/2)), ...
                    circshift((1:n)-n/2,round(n/2)) );
V = sum(sum( (K1.^2 + K2.^2) .* abs(Bh) ));
