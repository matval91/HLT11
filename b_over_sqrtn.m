%% plot ne
[shots, data]=HLT11_shots();
% figure();
% hold on;
for i=1:length(shots)
%     disp(shots(i))
     nel=gdat(shots(i), 'nel'); 
     b0=gdat_tcv(shots(i), 'b0', 'time_out', nel.t); 
     data(i).nel = nel; data(i).b0 =b0;
%     plot(nel.t, -1.*b0.data./sqrt(nel.data), 'DisplayName', num2str(shots(i)))
 end

for i=1:size(data,2)
    trange = data(i).t;
    if trange(1)==trange(2)
        continue
    end
    [~, ind_tmin] = min(abs(data(i).nel.t-trange(1))); data(i).ind_tmin=ind_tmin;
    [~, ind_tmax] = min(abs(data(i).nel.t-trange(2))); data(i).ind_tmax=ind_tmax;
    disp(data(i).nel.t(ind_tmin));     disp(data(i).nel.t(ind_tmax));
    b0ne = -1.*data(i).b0.data./sqrt(data(i).nel.data); avgb0ne=mean(b0ne(ind_tmin:ind_tmax));
    data(i).avgb0ne = avgb0ne;
end

figure(); hold on;
for i=1:size(data,2)
    data(i).ind_tmax=ind_tmax;data(i).ind_tmin=ind_tmin;
    try
        scatter(mean(data(i).nel.data(ind_tmin: ind_tmax)), data(i).avgb0ne*1e10);
    catch
        continue
    end
end