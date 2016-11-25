function x = secant_falsi(x0,x1,F,rtol,atol)
	fo = F(x0);
	for i=1:MAXIT
		fn = F(x1);
		s = fn*(x1-x0)/(fn-fo); % correction
		if (F(x1 - s)*fn < 0)
			x0 = x1;
			fo = fn;
		end
		x1 = x1 - s;
		if (abs(s) < max(atol,rtol*min(abs([x0;x1]))))
			x = x1;
			return;
		end
	end
end
