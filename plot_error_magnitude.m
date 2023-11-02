function plot_error_magnitude(s1,e1,e2,e3,noise_type)
figure
plot(s1-e1,'r','LineWidth',2)
hold on
plot(s1-e2,'g','LineWidth',2)
hold on
plot(s1-e3,'b','LineWidth',2)
title(['Plot of Error Magnitude vs Sample Index for ',noise_type],'FontSize',20,'FontWeight','bold')
xlabel('Samples (n)','FontSize',18,'FontWeight','bold')
ylabel('Error Magnitude','FontSize',18,'FontWeight','bold')
legend({'LMS','NLMS','RLS'},'FontSize',14)
end

