function Dist = Calculate_Dist_alpha1(DistType, Location1, Location2, f)

%pfp = (Location2' * f)';
Algeb = (Location2' * f)';
Algeb = Algeb .* Location1;
Dist = sum(Algeb, 1) .^ 2;

if strcmp(DistType, 'sampson')
    epl1 = f * Location1;
    epl2 = f' * Location2;
    Dist = Dist ./ (epl1(1,:).^2 + epl1(2,:).^2 + epl2(1,:).^2 + epl2(2,:).^2);
end