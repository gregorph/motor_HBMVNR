function plotEllipse(x, y, color)
    % Compute covariance
    C = cov(x.', y.');
    mu = mean([x.' y.'],1);
    [V, D] = eig(C);
    t = linspace(0, 2*pi, 100);
    % 95% CI scale factor
    chi2_val = sqrt(5.991); % chi2inv(0.95,2) ≈ 5.991
    a = (V*sqrt(D))*[cos(t); sin(t)]*chi2_val;
    plot(mu(1)+a(1,:), mu(2)+a(2,:), 'Color', color, 'LineWidth',1.5);
end