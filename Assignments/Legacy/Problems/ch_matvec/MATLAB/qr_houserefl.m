function Z = qr_houserefl(v)
% Use qr decomposition to find ONB of complement of span(v)
    [X,R] = qr(v);

    % Remove first column X(:,1) \in span(v)
    Z = X(:,2:end);
end
