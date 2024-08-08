clear R_indices L_indices ii ii_minus ii_plus iR iL box L_range R_range temp_map freq_map
clear temp_trace_v mode_limit MBO_modes_map min_mode max_mode m_range nu nu_pot var
% Build path to build map
L_indices=37;
R_indicies=37;
for ii=1:36
    ii_minus=37-ii;
    ii_plus=37+ii;
    if ii_plus<=72
        L_indices=[L_indices ii_plus ii_minus];
        R_indicies=[R_indicies ii_plus ii_minus];
    else
        L_indices=[L_indices ii_minus];
        R_indicies=[R_indicies ii_minus];
    end
end

freq_dist=15;
mode_limit=1; % Limit the +/- range of modes that can be compared to
MBO_modes_map=zeros(72,72,'double'); % Initialize mode map
MBO_modes_map(37,37)=46; % Set center 'known' value
% MBO_modes_map(:,:)=46; % Set center 'known' value
% MBO_modes_map(36,36)=46; % Set center 'known' value
% MBO_modes_map(36,35)=46; % Set center 'known' value

freq_map=zeros(length(L_indices),length(L_indices),'double');
freq_map(freq_map==0)=NaN;
% Loop over L indices
for iR=R_indicies
    % Loop over R indices
    for iL=L_indices
%         for ii=1:2
%             if ii==2,iL=i1;iR=i2;end
%             if ii==1,iL=i2;iR=i1;end
            if MBO_modes_map(iL,iR)==0
                % Create a range for the L indices
                L_range_mode=iL-mode_limit:iL+mode_limit;
                L_range_mode(L_range_mode<=0)=72-L_range_mode(L_range_mode<=0);
                L_range_mode(L_range_mode>72)=L_range_mode(L_range_mode>72)-72;
                L_range_nu=iL-1:iL+1;
                L_range_nu(L_range_nu<=0)=72-L_range_nu(L_range_nu<=0);
                L_range_nu(L_range_nu>72)=L_range_nu(L_range_nu>72)-72;
                % Create a range for the R indices
                R_range_mode=iR-mode_limit:iR+mode_limit;
                R_range_mode(R_range_mode<=0)=72-R_range_mode(R_range_mode<=0);
                R_range_mode(R_range_mode>72)=R_range_mode(R_range_mode>72)-72;
                R_range_nu=iR-1:iR+1;
                R_range_nu(R_range_nu<=0)=72-R_range_nu(R_range_nu<=0);
                R_range_nu(R_range_nu>72)=R_range_nu(R_range_nu>72)-72;
                % Create a temporary map from modes map within the box defined from
                % L and R indices
                temp_map=MBO_modes_map(L_range_mode,R_range_mode);
                temp_map(temp_map==0)=NaN; % Set any undefined values as NaNs
                % Define minimum mode allowed for this coordinate
                min_mode=max(temp_map(:))-mode_limit;
                % Define minimum mode allowed for this coordinate
                max_mode=min(temp_map(:))+mode_limit;
                % Define mode range for this coordinate
                if max_mode>min_mode,m_range=min_mode:max_mode;end
                if max_mode<min_mode,m_range=max_mode:min_mode;end
                if max_mode==min_mode,m_range=max_mode;end
                nu=freq_map(L_range_nu,R_range_nu);
                nu_max=max(nu(:));
                nu_min=min(nu(:));
                nu_pot=frequencies(m_range,iL,iR);
                nu_pot(abs(nu_pot-nu_max)>freq_dist)=[];
                nu_pot(abs(nu_pot-nu_min)>freq_dist)=[];
                [~,mode_inds]=ismember(nu_pot,frequencies(m_range,iL,iR));
                if ~isnan(nu_pot) % If there are good values to assign
                    mode_pot=m_range(mode_inds);
                    % Define a vector of scalar products within which to look for best mode
                    temp_trace_v=permute(Tr_vm(iL,iR,mode_pot)./Tr_v(iL,iR,mode_pot).*Tr_v_candidates(iL,iR,mode_pot),[3 1 2]);
                    [~,index]=max(temp_trace_v);
                    MBO_modes_map(iL,iR)=mode_pot(index);
                    freq_map(iL,iR)=frequencies(mode_pot(index),iL,iR);
                else % if no good options exist assign the average of the nearby points
                    freq_map(iL,iR)=sum(sum(nu(~isnan(nu))))/sum(sum(~isnan(nu)));
                end
            else % If mode is pre-assigned
                freq_map(iL,iR)=frequencies(MBO_modes_map(iL,iR),iL,iR);
            end
%         end
    end
end
freq_map(73,:)=freq_map(1,:);
freq_map(:,73)=freq_map(:,1);


figure;contourf(freq_map,50,'edgecolor','none');
colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');

freq_map9=repmat(freq_map,3,3);
freq_map9(freq_map9<0)=NaN;
good_data=~isnan(freq_map9);
[good_R_prop,good_L_prop]=meshgrid([[R_prop 180]-365 [R_prop 180] [R_prop 180]+365],[[R_prop 180]-365 [R_prop 180] [R_prop 180]+365]);
[full_axis_R,full_axis_L]=meshgrid(-180:179,-180:179);
freq_map=griddata(good_R_prop(good_data),good_L_prop(good_data),freq_map9(good_data),full_axis_R,full_axis_L);

freq_map=filter([0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1],1,freq_map);
axis_vec=-180:179;
figure;mesh(axis_vec(122:222),axis_vec(122:222),freq_map(122:222,122:222));
colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');
figure;contourf(axis_vec(122:222),axis_vec(122:222),freq_map(122:222,122:222),50,'edgecolor','none');
colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');

if true
figure;contourf(axis_vec-81,axis_vec(10:end)-81,freq_map(10:end,:),50,'edgecolor','none');
colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');
figure;mesh(axis_vec-81,axis_vec(10:end)-81,freq_map(10:end,:));
colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');
end


% freq_map(freq_map<0)=NaN;
% good_data=~isnan(freq_map);
% [good_R_prop,good_L_prop]=meshgrid([R_prop 180],[L_prop 180]);
% [full_axis_R,full_axis_L]=meshgrid(-180:179,-180:179);
% freq_map=griddata(good_R_prop(good_data),good_L_prop(good_data),freq_map(good_data),full_axis_R,full_axis_L);
% 
% freq_map=filter([0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1],1,freq_map);
% axis_vec=-180:179;
% figure;contourf(axis_vec(146:216),axis_vec(146:216),freq_map(146:216,146:216),50,'edgecolor','none');
% colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');
% 
% if false
% figure;contourf(axis_vec-81,axis_vec-81,freq_map,50,'edgecolor','none');
% colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');
% figure;contourf(axis_vec-81,axis_vec(10:end)-81,freq_map(10:end,:),50,'edgecolor','none');
% colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');
% end
