% initiate Tomlab environment
if (~exist('tomlablic', 'file'))
    if (strfind(lower(getenv('OS')), 'windows'))
        if (exist('C:\Program Files\tomlab\startup.m', 'file'))
            run('C:\Program Files\tomlab\startup');
        else
            error('Tomlab environment not found. Make sure to install Tomlab and run initTomlab or startup.');
        end
    else
        if (exist('/home/gbasler/tomlab/startup.m', 'file'))
            run('/home/gbasler/tomlab/startup');        
        elseif (exist('/package/devel/tomlab/startup.m', 'file'))
            run('/package/devel/tomlab/startup');
        elseif (exist('/usr/local/tomlab7.6/startup.m', 'file'))
            run('/usr/local/tomlab7.6/startup');
        elseif (exist('/usr/local/tomlab/startup.m', 'file'))
            run('/usr/local/tomlab/startup');
        else
            error('Tomlab environment not found. Make sure to install Tomlab and run initTomlab or startup.');
        end
    end
end
