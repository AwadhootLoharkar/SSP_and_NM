function [reflection, transmission, absorption] = RTA(energyStep, energy, damping, potentialPositions, potentialValues, delta, energyScale)
    % Calculate reflection, transmission, and absorption coefficients
    numE = length(energy);
    for kk = 1:numE
        [reflection(kk), transmission(kk), absorption(kk)] = RTA_iter(energyStep, energy(kk), damping, potentialPositions, potentialValues, delta, energyScale);
    end

    % Plot R(E), T(E), and A(E)
    figure;
    plot(energy, reflection, 'LineWidth', 2);
    hold on;
    plot(energy, transmission, 'LineWidth', 2);
    plot(energy, absorption, 'LineWidth', 2);
    xlabel('E (eV)', 'FontSize', 26);
    ylabel('R(E), T(E), A(E)', 'FontSize', 26);
    fontsize(gca, 22, 'points');
    ylim([0, 1]);
    set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
    set(gcf, 'Color', [0.7 0.7 0.7]);
    grid on;
    legend('R(E)', 'T(E)', 'A(E)', 'Location', 'Best');
end