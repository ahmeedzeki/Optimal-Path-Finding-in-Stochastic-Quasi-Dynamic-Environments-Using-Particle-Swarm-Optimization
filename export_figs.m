clear
clc
for NP = [2 4 6]
    for alpha = [0 100 10000]
        for npop = 50:50:500
            path_for_figs = ['./results_alpha_' num2str(alpha) '_np_' num2str(NP) '/npop_' num2str(npop) '.fig'];
            path_for_pdf = ['./results_alpha_' num2str(alpha) '_np_' num2str(NP) '/npop_' num2str(npop) '.pdf'];
            h = open(path_for_figs)
            xlabel('x')
            ylabel ('y')
            export_fig ('-transparent', '-pdf', path_for_pdf)
            close(h)
        end
    end
end
