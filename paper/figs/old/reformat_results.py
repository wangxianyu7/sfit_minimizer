
def format_event(event):
    #for line in event:
    #    print(line.strip())

    event_name = event[0]
    sfit = event[1]
    neldar = event[3]
    newton = event[4]
    bfgs = event[5]
    fits = [sfit, neldar, newton, bfgs]

    best_chi2 = None
    for fit in fits:
        chi2 = float(fit[32:43])
        if best_chi2 is None:
            best_chi2 = chi2
        elif chi2 < best_chi2:
            best_chi2 = chi2

    formatted_result = ''
    for fit in fits:
        formatted_result += '{0:<8} {1:12.2f} {2}\n'.format(event_name[0:8], best_chi2, fit.strip())

    return formatted_result


data = open('fits_results.txt', 'r')
outfile = open('fits_results_formatted.txt', 'w')

headstr = '{0:<9} {1:>11} '.format('Event', 'BestChi2')
headstr += '{0:15}'.format('Method')
headstr += '{0:8}'.format('Success')
headstr += '{0:>8}'.format('nfev')
headstr += '{0:>11}'.format('Chi2')
headstr += '{0:>11}'.format('t_0')
headstr += '{0:>10}'.format('u_0')
headstr += '{0:>9}'.format('t_E')
headstr += '{0:>9}'.format('njev')
outfile.write('{0}\n'.format(headstr))

event = None
for line in data.readlines():
    if event is None:
        event = [line]
    elif line[0:2] == 'KB':
        outfile.write(format_event(event))
        event = [line]
    else:
        event.append(line)

outfile.close()
data.close()