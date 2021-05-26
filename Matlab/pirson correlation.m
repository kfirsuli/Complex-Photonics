s1 = 25+1;
s2 = 35+1;
i1 = 4+1;
i2 = 10+1;

su = h(s1,i1)+h(s1,i2)+h(s2,i1)+h(s2,i2);
S_av = ((h(s1,i1)+h(s1,i2))*s1+(h(s2,i1)+h(s2,i2))*s2)/su;
I_va = ((h(s1,i1)+h(s2,i1))*i1+(h(s1,i2)+h(s2,i2))*i2)/su;

S_sdv = sqrt(((h(s1,i1)+h(s1,i2))*(s1-S_av)^2+(h(s2,i1)+h(s2,i2))*(s2-S_av)^2)/(su-1));
I_sdv = sqrt(((h(s1,i1)+h(s2,i1))*(i1-I_va)^2+(h(s1,i2)+h(s2,i2))*(i2-I_va)^2)/(su-1));

cor = h(s1,i1)*(s1-S_av)*(i1-I_va)+h(s1,i2)*(s1-S_av)*(i2-I_va)+h(s2,i1)*(s2-S_av)*(i1-I_va)+h(s2,i2)*(s2-S_av)*(i2-I_va);
cor = cor/S_sdv/I_sdv/(su-1)
%%
clear cor
for s1 = 15:45
    for s2 = 15:45
        for i1 = 1:15
            for i2 = 1:15
                su = h(s1,i1)+h(s1,i2)+h(s2,i1)+h(s2,i2);
S_av = ((h(s1,i1)+h(s1,i2))*s1+(h(s2,i1)+h(s2,i2))*s2)/su;
I_va = ((h(s1,i1)+h(s2,i1))*i1+(h(s1,i2)+h(s2,i2))*i2)/su;

S_sdv = sqrt(((h(s1,i1)+h(s1,i2))*(s1-S_av)^2+(h(s2,i1)+h(s2,i2))*(s2-S_av)^2)/(su-1));
I_sdv = sqrt(((h(s1,i1)+h(s2,i1))*(i1-I_va)^2+(h(s1,i2)+h(s2,i2))*(i2-I_va)^2)/(su-1));

cor(s1,s2,i1,i2) = h(s1,i1)*(s1-S_av)*(i1-I_va)+h(s1,i2)*(s1-S_av)*(i2-I_va)+h(s2,i1)*(s2-S_av)*(i1-I_va)+h(s2,i2)*(s2-S_av)*(i2-I_va);
cor(s1,s2,i1,i2) = cor(s1,s2,i1,i2)/S_sdv/I_sdv/(su-1);
            end
        end
    end
end
[C,I] = max(abs(cor(:)));
C
[s1,s2,i3,i4] = ind2sub(size(cor),I)
n = h(s1,i3)+h(s1,i4)+h(s2,i3)+h(s2,i4);
sqrt((1-C^2))/sqrt(h(s1,i3)+h(s1,i4)+h(s2,i3)+h(s2,i4)-2)


z = 0.5 * log((1+C)./(1-C));
zalpha = (-erfinv(0.35-1)) .* sqrt(2) ./ sqrt(n-3);
rlo = tanh(z-zalpha)-C
rup = tanh(z+zalpha)-C