function [greenFunction] = GreensFun(energyStep, position1, position2, energy, damping, energyScale)
    % Calculate the Green's function for given parameters and positions
    waveNumberOutside = sqrt((energy + 1i * damping) / energyScale);
    waveNumberInside = sqrt((energy + 1i * damping - energyStep) / energyScale);
    commonTerm = 1 / (2i * waveNumberInside);

    if position1 >= 0 && position2 >= 0
        greenFunction = commonTerm * (exp(1i * waveNumberInside * abs(position1 - position2)) + exp(1i * waveNumberInside * (position1 + position2)) * ((waveNumberInside - waveNumberOutside) / (waveNumberInside + waveNumberOutside)));
    elseif position1 < 0 && position2 >= 0 || position2 < 0 && position1 >= 0
        greenFunction = exp(-1i * waveNumberOutside * min(position1, position2) + 1i * waveNumberInside * max(position1, position2)) / (1i * (waveNumberOutside + waveNumberInside));
    end
end