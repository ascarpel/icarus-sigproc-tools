import uproot
import numpy as np
import sys

def plot_chanfft(plt, pdf, fftd, crate_label, event ):

    fig, ax = plt.subplots(1,2, figsize=(12,4) )

    for idx, plane_label in enumerate(['0', '1']):

        ax[idx].set_title(crate_label+"_plane"+plane_label+"_event"+str(event))
        out =ax[idx].imshow(fftd[plane_label], vmin=0, vmax=200, aspect=8)
        ax[idx].set_xlabel('Frequency (MHz)')
        ax[idx].set_ylabel('Ch. number')

        fig.colorbar(out, ax=ax[idx], label='Power', pad=0.02)

    pdf.savefig( fig  )
    plt.close(fig)

# ------------------------------------------------------------------------------

def plot_correlation(plt, pdf, corr, crate_label, event ):

    fig, ax = plt.subplots(1,2, figsize=(12,4) )

    for idx, plane_label in enumerate(['0', '1']):

        ax[idx].set_title(crate_label+"_plane"+plane_label+"_event"+str(event))
        out =ax[idx].imshow(corr[plane_label], vmin=-1, vmax=1, )
        ax[idx].set_xlabel('Ch number')
        ax[idx].set_ylabel('Ch. number')

        fig.colorbar(out, ax=ax[idx], label='Correlation', pad=0.02)

    pdf.savefig( fig  )
    plt.close(fig)

# ------------------------------------------------------------------------------

def plot_evd(plt, pdf, crate, crate_label, event ):

    fig, ax = plt.subplots(1,2, figsize=(12,4) )

    for idx, plane_label in enumerate(['0', '1']):

        ax[idx].set_title(crate_label+"_plane"+plane_label+"_event"+str(event))
        out =ax[idx].imshow(crate[plane_label], vmin=-4, vmax=4, aspect = 12)
        ax[idx].set_xlabel('Time ticks')
        ax[idx].set_ylabel('Ch. number')

        fig.colorbar(out, ax=ax[idx], label='ADC counts', pad=0.02)

    pdf.savefig( fig  )
    plt.close(fig)

# //////////////////////////////////////////////////////////////////////////////

def main():
    sigProcPath = "/Users/scarpell/Desktop/ICARUS/src/signal_processing/icarus-sigproc-tools"
    sys.path.insert(0,sigProcPath)

    PATHNAME       = "/Users/scarpell/Desktop/ICARUS/src/signal_processing/"
    PATHTODATA     = "data/noise/"
    RECOFILENAME   = PATHNAME + PATHTODATA + "/data_dl1_run822_201_20200109T010140_dl3_20200128T175402-DUMMY.root"

    # Below should be standard for the test data files currently available
    RECOFOLDERNAME = "Events"
    DAQNAME        = "raw::RawDigits_daq__DUMMY."

    ## Acquire the raw event #######################################################

    nchannels_minicrate = 576
    nchannels = 4608
    nchannels_minicrate = 576
    nboard = 64
    nchables = 32

    event = 1

    from sigproc_tools.sigproc_objects.rawdigit import RawDigit

    print("Opening file: ",RECOFILENAME)
    data_file = uproot.open(RECOFILENAME)

    print("Opening the folder contianing the RawDigits information: ",RECOFOLDERNAME)
    events_folder = data_file[RECOFOLDERNAME]

    rawdigits = RawDigit(events_folder,DAQNAME)
    rawWaveform = rawdigits.getWaveforms(event)

    ## Remove median baseline and RMS ##########################################

    from sigproc_tools.sigproc_functions.noiseProcessing import getPedestalsAndRMS
    waveforms, pedestal, rms = getPedestalsAndRMS(rawWaveform)

    # Study each of the 8 minicrate independently ###############################

    from sigproc_tools.sigproc_objects.minicrate import MiniCrate

    labels = [ 'crate'+str(i) for i in range(8) ]
    labels.remove('crate6') #not working

    # Save the crates in a dictionaly
    crates = { label : MiniCrate( int(label[-1]), waveforms ) for label in labels }

    # Now we print the relevant plots ##########################################

    from matplotlib.backends import backend_pdf
    import matplotlib.pyplot as plt

    pdf = backend_pdf.PdfPages("crates.pdf")

    for label in labels:

        print("Processing crate: " + label)
        crate = crates[label]

        plot_evd(plt, pdf,crate.plane , label, event )
        plot_correlation(plt, pdf, crate.get_corr() , label, event )
        fftd=crate.get_chanfft(1.25)
        plot_chanfft(plt, pdf, fftd , label, event )

        rmsd=crate.get_chanrms()

        # Now remove the coherent noise and replot everything ------------------
        crate.remove_coherent_noise(32)

        plot_evd(plt, pdf,crate.plane , label, event )
        plot_correlation(plt, pdf, crate.get_corr() , label, event )
        fftd_cnr=crate.get_chanfft(1.25)
        plot_chanfft(plt, pdf, fftd_cnr , label, event )

        rmsd_cnr=crate.get_chanrms()

        # Plot the comparison between rms --------------------------------------

        fig, ax = plt.subplots(1,2, figsize=(12,4) )

        for idx, plane_label in enumerate(['0', '1']):

            ax[idx].set_title(label+"_plane"+plane_label+"_event"+str(event))

            ax[idx].plot(np.linspace(0.0, crate.nchannels_minicrate/2, crate.nchannels_minicrate/2 ), rmsd[plane_label], label='Before CNR' )
            ax[idx].plot(np.linspace(0.0, crate.nchannels_minicrate/2, crate.nchannels_minicrate/2 ), rmsd_cnr[plane_label], label='After CNR' )
            #ax[0].set_yscale('Log')
            #ax[idx].set_xlim((0.0, 0.4))
            #ax[idx].set_ylim((0.1, 250))
            ax[idx].set_xlabel('Channel')
            ax[idx].set_ylabel('RMS (ADC)')
            ax[idx].legend()

        pdf.savefig( fig  )
        plt.close(fig)

        # Plot the comparison between ffts -------------------------------------

        fig, ax = plt.subplots(1,2, figsize=(12,4) )

        for idx, plane_label in enumerate(['0', '1']):

            ax[idx].set_title(label+"_plane"+plane_label+"_event"+str(event))

            ax[idx].plot(np.linspace(0.0, 1.25, fftd[plane_label].shape[1] ), np.mean(fftd[plane_label], axis=0), label='Noise' )
            ax[idx].plot(np.linspace(0.0, 1.25, fftd_cnr[plane_label].shape[1] ), np.mean(fftd_cnr[plane_label], axis=0), label='Denoised' )
            ax[idx].plot(np.linspace(0.0, 1.25, fftd_cnr[plane_label].shape[1] ), np.mean(fftd[plane_label], axis=0) - np.mean(fftd_cnr[plane_label], axis=0)  , label='Difference' )
            #ax[0].set_yscale('Log')
            ax[idx].set_xlim((0.0, 0.4))
            ax[idx].set_ylim((0.1, 250))
            ax[idx].set_xlabel('Frequency (MHz)')
            ax[idx].set_ylabel('Power')
            ax[idx].legend()

        pdf.savefig( fig  )
        plt.close(fig)

    pdf.close()

if __name__=='__main__':
    main()
