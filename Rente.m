% monatlicher Anfangsbeitragssatz
bmstart = 85;
% Gesambeitrag, Startwert
btstart = 0;
% Laufzeit in Jahren
lz = 38;
% Jaehrliche Erhoehung des monatlichen Beitragssatzes
je = 0.03;
% Tagesgeldkontozins der DKB
tknz = 0.02;
% Vermoegenssteuerzinssatz
vsz = 0.30;
z = 12;
tkez = (1+tknz/z)^z-1;

bm = bmstart;
bt = btstart;
tkbtn = btstart;
tkbtnw = btstart;
for yy = 1:lz
    % Gesamtbeitrag erhoeht sich jedes Jahr um 12*Monatsbeitrag
    bt = bt + 12*bm;
    % Tagesgeldkontoguthaben
    % monatlicher Zinszuwachs des Tagesgeldkontos, effectiv und nominal
    tkbtn0 = tkbtn;
    tkbtnw0 = tkbtnw;
    for mm = 1:12
        tkbtn0 = tkbtn0 + bmstart;
        tkbtn = (tkbtn + bmstart)*(1+tknz/12);
        tkbtnw0 = tkbtnw0 + bm;
        tkbtnw = (tkbtnw + bm)*(1+tknz/12);
    end
    dtkbtn = tkbtn-tkbtn0;
    dtkbtnw = tkbtnw-tkbtnw0;
    if dtkbtn > 800
        tkbtn = tkbtn - (dtkbtn-800)*vsz;
        tkbtnw = tkbtnw - (dtkbtnw-800)*vsz;
    end
    %fprintf('Jahr: %2u, monatl. Beitrag: %.2f\n',yy,bm)
    fprintf('Eingezahlt: %5.0f, mit Wachstum: %5.0f, Tagesgeldkontoguthaben: %8.2f, mit Wachstum: %8.2f, %8.2f (Zinsen)\n',yy*12*bmstart,bt,tkbtn,tkbtnw,dtkbtn)
    if yy ~= lz
        % Monatsbeitrag erhoeht sich nach jedem Jahr um 5%
	if mod(yy,3)==0
	   bm = bm*(1 + je);
	end
    end
end;

fprintf('\nBerechnung des eingezahlten Kapitals zur Rentenversicherung\n')
fprintf('Laufzeit: %u\nJaehrliche Erhoehung des monatlichen Beitragssatzes: %.3f%%\nAnfangsbeitragssatz: %6.2f\nBeitragssatz nach Ende der Laufzeit: %6.2f\n',lz,je,bmstart,bm)
fprintf('Gesamtbeitraege ohne jaehrliche Erhoehung: %.2f\n',lz*12*bmstart)
fprintf('Gesamtbeitraege mit jaehrliche Erhoehung:  %.2f\n',bt)

fprintf('\nAnlage auf Tagesgeldkonto\n')
fprintf('Nominaler Jahreszinssatz: %.5f\n',tknz)
fprintf('Effektiver Jahreszinssatz: %.5f\n',tkez)
fprintf('Gesamtbeitraege: %.2f\n',lz*12*bmstart)
fprintf('Gesamtguthaben bei Anlage auf Tagesgeldkonto, nominal:  %.2f\n',tkbtn)
