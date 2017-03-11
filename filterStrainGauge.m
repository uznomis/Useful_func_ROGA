function [test] = filterStrainGauge(test,cardNo,ch)
for i = 1:cardNo
    for j = 1:length(ch(i,:))
        if ch(i,j) == 0
            continue
        end
        test{i}(:,ch(i,j)) = ifilter(1:length(test{i}(:,ch(i,j))),...
            test{i}(:,ch(i,j)),...
    1.9565e-05,5.8695e-06,7,2,'Band-reject (notch)');
    end
end