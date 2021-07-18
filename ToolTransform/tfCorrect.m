function scan = tfCorrect(scan,tf)
% transfrom scan with the tf from base_laser to base_link

t = tf.translation;
if size(t,1) < size(t,2)
    t = t';
end
[y,~,~] = quat2angle(tf.rotation);
R = theta2R(y);

for i = 1:size(scan,2)
    if ~isempty(scan{i})
    scan{i}(2:3,:) = R*scan{i}(2:3,:) + t(1:2);
    end
end

end