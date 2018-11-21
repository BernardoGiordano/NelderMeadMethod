function res = Bz(X, z)
    I = X(1);
    R = X(2);
    Z = X(3);
    res = (I*R.^2)./(((Z-z).^2+R.^2).^(3/2)) + (I*R.^2)./(((Z+z).^2+R.^2).^(3/2));
end