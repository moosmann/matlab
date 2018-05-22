ca
% Create a figure and axes
vol = phantom3d( 200,200,200);
h = implay( vol );
f = h.Parent;
f.Visible = 'off';
ax = h.Visual.Axes;

t = h.DataSource.getTimeOfDisplayData;
t0 = h.DataSource.getOriginTime;
te = h.DataSource.getEndTime;
dt = h.DataSource.getSampleTimes;

h.DataSource.Controls.CurrentFrame;

txt_pos = [ceil(0.9*ax.Position(1:2)) 100 20];

% Add a text uicontrol to label the slider.
txt = uicontrol(f, 'Style','text',...
    'Position',txt_pos,...
    'String', sprintf('POS:%u', h.DataSource.Controls.CurrentFrame));

% Make figure visble after adding all components
f.Visible = 'on';

%set(f, 'ButtonDownFcn',@rotaxoff);
%set(f, 'WindowButtonDownFcn',@rotaxoff);
set(f,'WindowButtonDownFcn','test_rotaxoff(f,h,txt_pos)')
set(f,'ButtonDownFcn','test_rotaxoff(f,h,txt_pos)')
%set(f,'KeyPressFcn','test_rotaxoff(f,h,txt_pos)') does update ui but frame
% is not stepping anymore


%addlistener(h.DataSource.Controls,'CurrentFrame','PostSet',@(src,evnt)disp('Color changed'))
listener(h.DataSource.Controls,'CurrentFrame','PostSet',@(src,evnt)disp('Color changed'))
