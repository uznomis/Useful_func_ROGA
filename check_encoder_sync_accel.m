% check encoder alignment for Accel
figure;
for i = 1:4
hold on
plot(Accel(:,5*i-4),Accel(:,5*i-3));
hold off
end