function  Plot2D(truevalue, pf, marker)
pause(0.01);
plot(truevalue(:,1), truevalue(:,2),'.');
hold on;
plot(pf(:,1), pf(:,2), marker, 'MarkerFace', 'r');
% xlim([0 max(truevalue(:,1)) + 0.1]);
% ylim([0 max(truevalue(:,2)) + 0.1]);
xlim([0 max(pf(:,1)) + 0.1]);
ylim([0 max(pf(:,2)) + 0.1]);
xlabel('f_1');ylabel('f_2');
hold off;
end