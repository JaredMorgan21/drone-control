function zd = predict_trajectory(t, reading)    
    syms time
    z = reading(time);
    v = diff(z);
    a = diff(z,2);

    z = subs(z, time, t);
    v = subs(v, time, t);
    a = subs(a, time, t);

    zd = z + 0.1*v + 0.01*a;
    zd = zd(1:12);
    zd = round(zd, 2);
end