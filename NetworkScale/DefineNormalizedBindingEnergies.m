function DefineNormalizedBindingEnergies

global ka kd b D length_conversion...
    eaStar edStar

eaStar = -log(ka*(b*length_conversion)^2/D);
edStar = -log(kd*(b*length_conversion)^2/D);

end