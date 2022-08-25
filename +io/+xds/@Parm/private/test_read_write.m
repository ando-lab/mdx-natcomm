write_xparm(load_xparm('XPARM.XDS'),'test.txt')
[~,a] = unix('tee < XPARM.XDS');
[~,b] = unix('tee < test.txt');
if length(a)==length(b) && all(a==b)
    disp('read and write works as expected');
else
    disp('read and write did NOT work correctly...');
end