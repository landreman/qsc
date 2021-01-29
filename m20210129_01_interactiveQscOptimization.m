function m20210129_01_interactiveQscOptimization()

% Absolute path to the qsc executable:
qsc_executable = '/Users/mattland/qsc/bin/xqsc';

% Qsc input file to use as a template:
qsc_template = '/Users/mattland/Box Sync/work21/20210129-01-qsc_optimize_QH4/qsc_in.20210129-01-QH4_006_fromScan';

% Input file to write:
qsc_input = '/Users/mattland/Box Sync/work21/20210129-01-qsc_optimize_QH4/qsc_in.interactive';

template = splitlines(fileread(qsc_template));

weight_XY2 = 0.0;
weight_XY2Prime = 0.0;
weight_XY3 = 0.0;
weight_XY3Prime = 0.0;
weight_B20 = 0.0;
weight_grad_grad_B = 0.0;
weight_r_singularity = 0.0;
weight_iota = 0.0;
weight_well = 0.0;
weight_R0 = 0.0;

target_iota = 0.0;
target_well = 0.0;
target_R0 = 0.0;

first_update = true;
number_of_calls = 0;

B20_term_handle = 0;
XY2_term_handle = 0;
XY2Prime_term_handle = 0;
XY3_term_handle = 0;
XY3Prime_term_handle = 0;
grad_grad_B_term_handle = 0;
r_singularity_term_handle = 0;
iota_term_handle = 0;
well_term_handle = 0;
R0_term_handle = 0;
objective_function_handle = 0;

eta_bar_handle = 0;
B2c_handle = 0;
R0c_handle = 0;
Z0s_handle = 0;
iota_handle = 0;
B20_handle = 0;
r_singularity_handle = 0;
L_grad_B_handle = 0;
L_grad_grad_B_handle = 0;
well_handle = 0;

f = figure('Visible','off','Units','pixels','Position',[0,0,1600,800]);

%label_margin=25;
label_margin=20;
margin = 44;
height=665;
label_left = 10;
label_width=200;
text_left = 110;
text_width=100;
button_width = 200;
button_height = 25;
button_margin = 17;
slider_left = 15;
slider_width = 200;
slider_height = 20;
double_slider_margin = 34;
button_fontsize = 10;
eventName = 'PostSet';
weight_min = -3;
weight_max = 3;

B20_button = uicontrol('Style','checkbox','Units','pixels','value',false,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include B20','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_B20 = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_B20,'Value',eventName,@slider_B20_callback);
label_B20 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_B20','horizontalalignment','left');
text_B20 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_B20_callback);

height = height - margin;
XY2_button = uicontrol('Style','checkbox','Units','pixels','value',true,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include XY2','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_XY2 = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_XY2,'Value',eventName,@slider_XY2_callback);
label_XY2 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_XY2','horizontalalignment','left');
text_XY2 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_XY2_callback);

height = height - margin;
XY2Prime_button = uicontrol('Style','checkbox','Units','pixels','value',true,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include XY2Prime','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_XY2Prime = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_XY2Prime,'Value',eventName,@slider_XY2Prime_callback);
label_XY2Prime = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_XY2Prime','horizontalalignment','left');
text_XY2Prime = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_XY2Prime_callback);

height = height - margin;
XY3_button = uicontrol('Style','checkbox','Units','pixels','value',false,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include XY3','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_XY3 = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_XY3,'Value',eventName,@slider_XY3_callback);
label_XY3 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_XY3','horizontalalignment','left');
text_XY3 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_XY3_callback);

height = height - margin;
XY3Prime_button = uicontrol('Style','checkbox','Units','pixels','value',false,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include XY3Prime','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_XY3Prime = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_XY3Prime,'Value',eventName,@slider_XY3Prime_callback);
label_XY3Prime = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_XY3Prime','horizontalalignment','left');
text_XY3Prime = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_XY3Prime_callback);

height = height - margin;
grad_grad_B_button = uicontrol('Style','checkbox','Units','pixels','value',false,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include grad_grad_B','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_grad_grad_B = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_grad_grad_B,'Value',eventName,@slider_grad_grad_B_callback);
label_grad_grad_B = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_grad_grad_B','horizontalalignment','left');
text_grad_grad_B = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_grad_grad_B_callback);

height = height - margin;
r_singularity_button = uicontrol('Style','checkbox','Units','pixels','value',false,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include r_singularity','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_r_singularity = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_r_singularity,'Value',eventName,@slider_r_singularity_callback);
label_r_singularity = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_r_singularity','horizontalalignment','left');
text_r_singularity = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_r_singularity_callback);

% Beginning of terms that involve a target value:

height = height - margin;
iota_button = uicontrol('Style','checkbox','Units','pixels','value',false,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include iota','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_iota = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_iota,'Value',eventName,@slider_iota_callback);
label_iota = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_iota','horizontalalignment','left');
text_iota = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_iota_callback);
height = height - double_slider_margin;
slider_iota_target = uicontrol('Style','slider','Min',-2,'Max',2,'Value',-0.4,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_iota_target,'Value',eventName,@slider_iota_target_callback);
label_iota_target = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','target_iota','horizontalalignment','left');
text_iota_target = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_iota_target_callback);

height = height - margin;
well_button = uicontrol('Style','checkbox','Units','pixels','value',false,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include well','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_well = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_well,'Value',eventName,@slider_well_callback);
label_well = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_well','horizontalalignment','left');
text_well = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_well_callback);
height = height - double_slider_margin;
slider_well_target = uicontrol('Style','slider','Min',-150,'Max',0,'Value',0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_well_target,'Value',eventName,@slider_well_target_callback);
label_well_target = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','target_well','horizontalalignment','left');
text_well_target = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_well_target_callback);

height = height - margin;
R0_button = uicontrol('Style','checkbox','Units','pixels','value',false,'Position',[label_left,height+label_margin,button_width,button_height],'String','Include R0','Callback',@button_callback,'fontsize',button_fontsize);
height = height - button_margin;
slider_R0 = uicontrol('Style','slider','Min',weight_min,'Max',weight_max,'Value',0.0,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_R0,'Value',eventName,@slider_R0_callback);
label_R0 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','weight_R0','horizontalalignment','left');
text_R0 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_R0_callback);
height = height - double_slider_margin;
slider_R0_target = uicontrol('Style','slider','Min',-2,'Max',2,'Value',-0.4,'Units','pixels','Position',[slider_left,height,slider_width,slider_height]);
addlistener(slider_R0_target,'Value',eventName,@slider_R0_target_callback);
label_R0_target = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','target_R0','horizontalalignment','left');
text_R0_target = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_R0_target_callback);

set_text_XY2()
set_text_XY2Prime()
set_text_XY3()
set_text_XY3Prime()
set_text_B20()
set_text_grad_grad_B()
set_text_r_singularity()

set_text_iota()
set_text_iota_target()
set_text_well()
set_text_well_target()
set_text_R0()
set_text_R0_target()

update()

set(f,'Visible','on')

    function button_callback(source,callbackdata)
        update()
    end

    function slider_XY2_callback(source,callbackdata)
        weight_XY2 = slider_XY2.Value;
        set_text_XY2()
        update()
    end
    function slider_XY2Prime_callback(source,callbackdata)
        weight_XY2Prime = slider_XY2Prime.Value;
        set_text_XY2Prime()
        update()
    end
    function slider_XY3_callback(source,callbackdata)
        weight_XY3 = slider_XY3.Value;
        set_text_XY3()
        update()
    end
    function slider_XY3Prime_callback(source,callbackdata)
        weight_XY3Prime = slider_XY3Prime.Value;
        set_text_XY3Prime()
        update()
    end
    function slider_B20_callback(source,callbackdata)
        weight_B20 = slider_B20.Value;
        set_text_B20()
        update()
    end
    function slider_grad_grad_B_callback(source,callbackdata)
        weight_grad_grad_B = slider_grad_grad_B.Value;
        set_text_grad_grad_B()
        update()
    end
    function slider_r_singularity_callback(source,callbackdata)
        weight_r_singularity = slider_r_singularity.Value;
        set_text_r_singularity()
        update()
    end
    function slider_iota_callback(source,callbackdata)
        weight_iota = slider_iota.Value;
        set_text_iota()
        update()
    end
    function slider_iota_target_callback(source,callbackdata)
        target_iota = slider_iota_target.Value;
        set_text_iota_target()
        update()
    end
    function slider_well_callback(source,callbackdata)
        weight_well = slider_well.Value;
        set_text_well()
        update()
    end
    function slider_well_target_callback(source,callbackdata)
        target_well = slider_well_target.Value;
        set_text_well_target()
        update()
    end
    function slider_R0_callback(source,callbackdata)
        weight_R0 = slider_R0.Value;
        set_text_R0()
        update()
    end
    function slider_R0_target_callback(source,callbackdata)
        target_R0 = slider_R0_target.Value;
        set_text_R0_target()
        update()
    end

    function text_XY2_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_XY2.Max]),slider_XY2.Min]);
            weight_XY2 = val;
            slider_XY2.Value = weight_XY2;
            % The above line does not cause the slider callback to be called.
            set_text_XY2()
            update()
        else
            set_text_XY2()
        end
    end

    function text_XY2Prime_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_XY2Prime.Max]),slider_XY2Prime.Min]);
            weight_XY2Prime = val;
            slider_XY2Prime.Value = weight_XY2Prime;
            % The above line does not cause the slider callback to be called.
            set_text_XY2Prime()
            update()
        else
            set_text_XY2Prime()
        end
    end

    function text_XY3_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_XY3.Max]),slider_XY3.Min]);
            weight_XY3 = val;
            slider_XY3.Value = weight_XY3;
            % The above line does not cause the slider callback to be called.
            set_text_XY3()
            update()
        else
            set_text_XY3()
        end
    end

    function text_XY3Prime_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_XY3Prime.Max]),slider_XY3Prime.Min]);
            weight_XY3Prime = val;
            slider_XY3Prime.Value = weight_XY3Prime;
            % The above line does not cause the slider callback to be called.
            set_text_XY3Prime()
            update()
        else
            set_text_XY3Prime()
        end
    end

    function text_B20_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_B20.Max]),slider_B20.Min]);
            weight_B20 = val;
            slider_B20.Value = weight_B20;
            % The above line does not cause the slider callback to be called.
            set_text_B20()
            update()
        else
            set_text_B20()
        end
    end

    function text_grad_grad_B_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_grad_grad_B.Max]),slider_grad_grad_B.Min]);
            weight_grad_grad_B = val;
            slider_grad_grad_B.Value = weight_grad_grad_B;
            % The above line does not cause the slider callback to be called.
            set_text_grad_grad_B()
            update()
        else
            set_text_grad_grad_B()
        end
    end

    function text_r_singularity_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_r_singularity.Max]),slider_r_singularity.Min]);
            weight_r_singularity = val;
            slider_r_singularity.Value = weight_r_singularity;
            % The above line does not cause the slider callback to be called.
            set_text_r_singularity()
            update()
        else
            set_text_r_singularity()
        end
    end

    function text_iota_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_iota.Max]),slider_iota.Min]);
            weight_iota = val;
            slider_iota.Value = weight_iota;
            % The above line does not cause the slider callback to be called.
            set_text_iota()
            update()
        else
            set_text_iota()
        end
    end

    function text_iota_target_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_iota_target.Max]),slider_iota_target.Min]);
            target_iota = val;
            slider_iota_target.Value = target_iota;
            % The above line does not cause the slider callback to be called.
            set_text_iota_target()
            update()
        else
            set_text_iota_target()
        end
    end

    function text_well_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_well.Max]),slider_well.Min]);
            weight_well = val;
            slider_well.Value = weight_well;
            % The above line does not cause the slider callback to be called.
            set_text_well()
            update()
        else
            set_text_well()
        end
    end

    function text_well_target_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_well_target.Max]),slider_well_target.Min]);
            target_well = val;
            slider_well_target.Value = target_well;
            % The above line does not cause the slider callback to be called.
            set_text_well_target()
            update()
        else
            set_text_well_target()
        end
    end

    function text_R0_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_R0.Max]),slider_R0.Min]);
            weight_R0 = val;
            slider_R0.Value = weight_R0;
            % The above line does not cause the slider callback to be called.
            set_text_R0()
            update()
        else
            set_text_R0()
        end
    end

    function text_R0_target_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_R0_target.Max]),slider_R0_target.Min]);
            target_R0 = val;
            slider_R0_target.Value = target_R0;
            % The above line does not cause the slider callback to be called.
            set_text_R0_target()
            update()
        else
            set_text_R0_target()
        end
    end

    function set_text_XY2
        text_XY2.String = sprintf('%g',weight_XY2);
    end
    function set_text_XY2Prime
        text_XY2Prime.String = sprintf('%g',weight_XY2Prime);
    end
    function set_text_XY3
        text_XY3.String = sprintf('%g',weight_XY3);
    end
    function set_text_XY3Prime
        text_XY3Prime.String = sprintf('%g',weight_XY3Prime);
    end
    function set_text_B20
        text_B20.String = sprintf('%g',weight_B20);
    end
    function set_text_grad_grad_B
        text_grad_grad_B.String = sprintf('%g',weight_grad_grad_B);
    end
    function set_text_r_singularity
        text_r_singularity.String = sprintf('%g',weight_r_singularity);
    end

    function set_text_iota
        text_iota.String = sprintf('%g',weight_iota);
    end
    function set_text_iota_target
        text_iota_target.String = sprintf('%g',target_iota);
    end

    function set_text_well
        text_well.String = sprintf('%g',weight_well);
    end
    function set_text_well_target
        text_well_target.String = sprintf('%g',target_well);
    end

    function set_text_R0
        text_R0.String = sprintf('%g',weight_R0);
    end
    function set_text_R0_target
        text_R0_target.String = sprintf('%g',target_R0);
    end


    function update()
        fprintf('Writing file\n');
        % Write qsc input file:
        fid = fopen(qsc_input,'w');
        if fid < 0
            error(['Unable to open file ',input_filename,' for writing.'])
        end
        for j = 1:numel(template)
            line = template{j};
            %fprintf('Handling line: %s\n',line);
            % Remove any lines containing "weight":
            %if (~contains(line, 'weight')) && (~isempty(line))
            if (~contains(line, 'weight') ...
                    && (~contains(line, 'target_iota')) ...
                    && (~contains(line, 'max_d2_volume_d_psi2')) ...
                    && (~contains(line, 'min_R0')))
                fprintf(fid, [line, '\n']);
            end
        end
        
        fprintf(fid,'target_iota = %g\n', target_iota);
        fprintf(fid,'max_d2_volume_d_psi2 = %g\n', target_well);
        fprintf(fid,'min_R0 = %g\n', target_R0);
        
        val = -1;
        if (XY2_button.Value == XY2_button.Max)
            val = 10 ^ weight_XY2;
        end
        fprintf(fid,'weight_XY2 = %g\n', val);
        
        val = -1;
        if (XY2Prime_button.Value == XY2Prime_button.Max)
            val = 10 ^ weight_XY2Prime;
        end
        fprintf(fid,'weight_XY2Prime = %g\n', val);

        val = -1;
        if (XY3_button.Value == XY3_button.Max)
            val = 10 ^ weight_XY3;
        end
        fprintf(fid,'weight_XY3 = %g\n', val);
        
        val = -1;
        if (XY3Prime_button.Value == XY3Prime_button.Max)
            val = 10 ^ weight_XY3Prime;
        end
        fprintf(fid,'weight_XY3Prime = %g\n', val);

        val = -1;
        if (B20_button.Value == B20_button.Max)
            val = 10 ^ weight_B20;
        end
        fprintf(fid,'weight_B20 = %g\n', val);

        val = -1;
        if (grad_grad_B_button.Value == grad_grad_B_button.Max)
            val = 10 ^ weight_grad_grad_B;
        end
        fprintf(fid,'weight_grad_grad_B = %g\n', val);

        val = -1;
        if (r_singularity_button.Value == r_singularity_button.Max)
            val = 10 ^ weight_r_singularity;
        end
        fprintf(fid,'weight_r_singularity = %g\n', val);

        val = -1;
        if (iota_button.Value == iota_button.Max)
            val = 10 ^ weight_iota;
        end
        fprintf(fid,'weight_iota = %g\n', val);

        val = -1;
        if (well_button.Value == well_button.Max)
            val = 10 ^ weight_well;
        end
        fprintf(fid,'weight_d2_volume_d_psi2 = %g\n', val);

        val = -1;
        if (R0_button.Value == R0_button.Max)
            val = 10 ^ weight_R0;
        end
        fprintf(fid,'weight_R0 = %g\n', val);

        fclose(fid);
        
        % Run the executable:
        fprintf([qsc_executable, ' ', qsc_input, '\n']);
        system([qsc_executable, ' "', qsc_input, '"']);
        number_of_calls = number_of_calls + 1;
        
        output_filename = [replace(qsc_input, 'qsc_in', 'qsc_out'), '.nc'];
        
        eta_bar = ncread(output_filename, 'iter_eta_bar');
        B2c = ncread(output_filename, 'iter_B2c');
        R0c = ncread(output_filename, 'iter_R0c');
        Z0s = ncread(output_filename, 'iter_Z0s');
        
        iota = ncread(output_filename, 'iter_iota');
        B20_variation = ncread(output_filename, 'iter_B20_variation');
        r_singularity = ncread(output_filename, 'iter_r_singularity');
        d2_volume_d_psi2 = ncread(output_filename, 'iter_d2_volume_d_psi2');
        L_grad_B = ncread(output_filename, 'iter_min_L_grad_B');
        L_grad_grad_B = ncread(output_filename, 'iter_min_L_grad_grad_B');
        elongation = ncread(output_filename, 'iter_max_elongation');
        
        objective_function = ncread(output_filename, 'iter_objective_function');
        B20_term = ncread(output_filename, 'iter_B20_term');
        iota_term = ncread(output_filename, 'iter_iota_term');
        R0_term = ncread(output_filename, 'iter_R0_term');
        well_term = ncread(output_filename, 'iter_d2_volume_d_psi2_term');
        XY2_term = ncread(output_filename, 'iter_XY2_term');
        XY2Prime_term = ncread(output_filename, 'iter_XY2Prime_term');
        XY3_term = ncread(output_filename, 'iter_XY3_term');
        XY3Prime_term = ncread(output_filename, 'iter_XY3Prime_term');
        grad_grad_B_term = ncread(output_filename, 'iter_grad_grad_B_term');
        r_singularity_term = ncread(output_filename, 'iter_r_singularity_term');
        
        nfourier = size(R0c, 1);

        %{
        nrows = 2;
        ncols = 2;
        
        subplot(nrows, ncols, 1)
        semilogy(objective_function, '.-');
        
        subplot(nrows, ncols, 2)
        semilogy(iota, '.-');
        
        subplot(nrows, ncols, 3)
        semilogy(r_singularity, '.-');
        %}
        
        if first_update
            first_update = false;
            
            plot_width = 220;
            plot_height = 185;
            plot_left = 250;
            plot_bottom = 35;
            plot_horiz_space = 265;
            plot_vert_space = 235;
            markersize = 4;
            
            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 0, plot_bottom + plot_vert_space * 2, plot_width, plot_height]);
            B20_term_handle = semilogy(B20_term, '.-', 'displayname', 'B20 term');
            hold on
            XY2_term_handle = semilogy(XY2_term, '.-', 'displayname', 'XY2 term');
            XY2Prime_term_handle = semilogy(XY2Prime_term, '.-', 'displayname', 'XY2Prime term');
            XY3_term_handle = semilogy(XY3_term, '.-', 'displayname', 'XY3 term');
            XY3Prime_term_handle = semilogy(XY3Prime_term, '.-', 'displayname', 'XY3Prime term');
            grad_grad_B_term_handle = semilogy(grad_grad_B_term, '.-', 'displayname', 'grad grad B term');
            r_singularity_term_handle = semilogy(r_singularity_term, '.-', 'displayname', 'r singularity term');
            iota_term_handle = semilogy(iota_term, 's-', 'markersize', markersize, 'displayname', 'iota term');
            well_term_handle = semilogy(well_term, 's-', 'markersize', markersize, 'displayname', 'well term');
            R0_term_handle = semilogy(R0_term, 's-', 'markersize', markersize, 'displayname', 'R0 term');
            objective_function_handle = semilogy(objective_function, ':k', 'displayname','total objective');
            ylim([1e-6, 1e4])
            legend('position',[0.3740    0.7717    0.0859    0.1815])
            
            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 2, plot_bottom + plot_vert_space * 2, plot_width, plot_height]);
            eta_bar_handle = plot(eta_bar, '.-');
            title('eta bar')
            ylim([-2.5, 2.5])

            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 3, plot_bottom + plot_vert_space * 2, plot_width, plot_height]);
            B2c_handle = plot(B2c, '.-');
            title('B2c')
            ylim([-3, 3])

            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 0, plot_bottom + plot_vert_space * 1, plot_width, plot_height]);
            R0c_handle = semilogy(abs(R0c'), '.-');
            title('R0c')
            ylim([1e-6, 1.2])

            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 1, plot_bottom + plot_vert_space * 1, plot_width, plot_height]);
            Z0s_handle = semilogy(abs(Z0s'), '.-');
            title('Z0s')
            ylim([1e-6, 1.2])

            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 2, plot_bottom + plot_vert_space * 1, plot_width, plot_height]);
            iota_handle = plot(iota, '.-');
            title('iota')
            ylim([-2, 2])

            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 3, plot_bottom + plot_vert_space * 1, plot_width, plot_height]);
            B20_handle = plot(B20_variation, '.-');
            title('B20 variation')
            ylim([0, 4])

            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 0, plot_bottom + plot_vert_space * 0, plot_width, plot_height]);
            r_singularity_handle = plot(r_singularity, '.-');
            title('r singularity')
            ylim([0, 0.6])
            xlabel('iteration')

            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 1, plot_bottom + plot_vert_space * 0, plot_width, plot_height]);
            L_grad_B_handle = plot(L_grad_B, '.-', 'displayname', 'L grad B');
            hold on
            L_grad_grad_B_handle = plot(L_grad_grad_B, '.-', 'displayname', 'L grad grad B');
            ylim([0, 0.6])
            legend show
            xlabel('iteration')

            ax = axes('Units','pixels','Position',[plot_left + plot_horiz_space * 2, plot_bottom + plot_vert_space * 0, plot_width, plot_height]);
            well_handle = plot(d2_volume_d_psi2, '.-');
            title('V"')
            ylim([-300, 1500])
            xlabel('iteration')

        else
            set(B20_term_handle, 'YData', B20_term);
            set(XY2_term_handle, 'YData', XY2_term);
            set(XY2Prime_term_handle, 'YData', XY2Prime_term);
            set(XY3_term_handle, 'YData', XY3_term);
            set(XY3Prime_term_handle, 'YData', XY3Prime_term);
            set(grad_grad_B_term_handle, 'YData', grad_grad_B_term);
            set(r_singularity_term_handle, 'YData', r_singularity_term);
            set(iota_term_handle, 'YData', iota_term);
            set(well_term_handle, 'YData', well_term);
            set(R0_term_handle, 'YData', R0_term);
            set(objective_function_handle, 'YData', objective_function);

            set(eta_bar_handle, 'YData', eta_bar);
            set(B2c_handle, 'YData', B2c);
            for j = 1:nfourier
                set(R0c_handle(j), 'YData', abs(R0c(j,:)));
                set(Z0s_handle(j), 'YData', abs(Z0s(j,:)));
            end
            set(iota_handle, 'YData', iota);
            set(B20_handle, 'YData', B20_variation);
            set(r_singularity_handle, 'YData', r_singularity);
            set(L_grad_B_handle, 'YData', L_grad_B);
            set(L_grad_grad_B_handle, 'YData', L_grad_grad_B);
            set(well_handle, 'YData', d2_volume_d_psi2);
            
            refreshdata
        end
        
        
    end


end