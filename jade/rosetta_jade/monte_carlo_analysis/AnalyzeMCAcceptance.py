
from collections import defaultdict
from jade.basic import path as path_tools
import numpy
from jade.basic.plotting.MakeFigure import MakeFigure








### I think I'll be deprecating this old baby in favor of ipython notebooks and pandas.

class AnalyzeMCAcceptance:
    """
    Calculates Monte Carlo acceptance and plots data by reading PDB files.
    """
    def __init__(self, analysis_infos, root_dataset_dir="datasets/pdblists",data_outdir="data", plot_outdir="plots/mc_benchmarks"):

        #Main Analysis classes that store paths to what we need
        self.analysis_infos = analysis_infos

        #Main paths for outputing and inputting outside of the analysis classes
        self.root_dataset_dir = root_dataset_dir
        self.data_outdir = data_outdir
        self.plot_outdir = plot_outdir

        #All of our resultant data.
        self.all_data = defaultdict()

    def analyze(self):
        for bm_info in self.analysis_infos:
            bm_data = MCData(self.root_dataset_dir, bm_info)
            self.all_data[bm_info.get_exp()] = bm_data

        self.plot_cum_acceptance()
        self.plot_cycle_by_e()



    def plot_cycle_by_e(self):
        def add_data(maker, all_models, include_native, final_e, range_only = None):
                i = 0
                for cy_data in all_models:
                    i+=1
                    assert isinstance(cy_data, CycleData)


                    y = []

                    if range_only:
                        x = [z for z in range(range_only[0], range_only[1] + 1)]
                    elif include_native:
                        x = [z for z in range(-1, 152)]
                        y.append(cy_data.native_e)
                        y.append(cy_data.start_e)
                    else:
                        x = [z for z in range(0, 152)]
                        y.append(cy_data.start_e)



                    if range_only:
                        r = range(range_only[0], range_only[1] )
                        #print repr(r)
                    else:
                        r = cy_data.values.keys()

                    if final_e:
                        y.extend([cy_data.get_final_e(ii) for ii in r])
                    else:
                        y.extend([cy_data.get_mc_e(ii) for ii in r])

                    #print str(len(x))
                    #print str(len(y))

                    y.append(cy_data.get_protocol_final_e())

                    maker.add_data(x, y, "model_"+repr(i))

        for exp in self.all_data:
            root_exp_dir = path_tools.get_make_get_dirs(self.plot_outdir, [exp, "bm_pdbs"])
            mc_data = self.all_data[exp]
            assert isinstance(mc_data, MCData)


            for pdb_id in mc_data.data.keys():

                all_models = mc_data.get(pdb_id)
                make_figure = MakeFigure(2, 1)

                #root_pdbid_dir = path_tools.get_make_get_dirs(root_exp_dir, [pdb_id])


                #### Plot not including antigen ####

                ###### Pre MC ######

                add_data(make_figure, all_models, include_native=True, final_e=False)
                make_figure.fill_subplot("Pre MC E "+pdb_id.split("_")[0], make_figure.labels, y_axis_label= "REU" , add_legend=False)



                ###### Post MC ######
                add_data(make_figure, all_models, include_native=True, final_e=True)
                make_figure.fill_subplot("Post MC E "+pdb_id.split("_")[0], make_figure.labels, x_axis_label="Cycle", y_axis_label="REU", add_legend=False)




                make_figure.set_y_scale('symlog')
                make_figure.add_grid()
                make_figure.save_plot(root_exp_dir+"/"+pdb_id.split("_")[0]+".cycle-vs-mc_e_log.with_native.pdf")


                #### Plot including antigen ####
                make_figure = MakeFigure(2, 1)

                ###### Pre MC ######
                add_data(make_figure, all_models, include_native=False, final_e=False)

                make_figure.fill_subplot("Pre MC E "+pdb_id.split("_")[0], make_figure.labels, y_axis_label= "REU" , add_legend=False)



                ###### Post MC ######
                add_data(make_figure, all_models, include_native=False, final_e=True)

                #Make first subplot
                make_figure.fill_subplot("Post MC E "+pdb_id.split("_")[0], make_figure.labels, x_axis_label="Cycle", y_axis_label="REU", add_legend=False)
                make_figure.set_y_scale('symlog')
                make_figure.add_grid()
                make_figure.save_plot(root_exp_dir+"/"+pdb_id.split("_")[0] +".cycle-vs-mc_e_log.no_native.pdf")

                make_figure = MakeFigure(1, 1)
                add_data(make_figure, all_models, include_native=False, final_e=True, range_only=[50, 150])
                make_figure.fill_subplot("Post MC E "+ pdb_id.split("_")[0] , make_figure.labels, x_axis_label="Cycle", y_axis_label="REU", add_legend=False)
                make_figure.add_grid()
                make_figure.save_plot(root_exp_dir+"/"+pdb_id.split("_")[0] +".cycle-vs-mc_e_From50.pdf")


    def plot_cum_acceptance(self):
        def add_data(maker, all_models):
                i = 0
                for cy_data in all_models:
                    i+=1
                    assert isinstance(cy_data, CycleData)

                    y = []

                    x = [z for z in range(1, len(cy_data.values.keys()) +1)]
                    #print repr(x)
                    r = cy_data.values.keys()
                    y.extend([cy_data.get_cum_acceptance(ii) for ii in r])

                    maker.add_data(x, y, "model_"+repr(i))


        make_exp_figure = MakeFigure(1, 1)


        x_values = None
        for bm_info in self.analysis_infos:
            exp = bm_info.get_exp()

            make_all_pdbids_figure = MakeFigure(1, 1)
            root_exp_dir = path_tools.get_make_get_dirs(self.plot_outdir, [exp, "bm_pdbs"])
            mc_data = self.all_data[exp]
            all_y_values = []

            for pdb_id in mc_data.data.keys():

                all_models = mc_data.get(pdb_id)
                make_figure = MakeFigure(1, 1)

                add_data(make_figure, all_models)
                make_figure.fill_subplot("Cumulative MC Acceptance "+pdb_id.split("_")[0], make_figure.labels, x_axis_label="Cycle", y_axis_label="Cummulative Accepts", add_legend=False)
                make_figure.add_grid()
                make_figure.save_plot(root_exp_dir+"/"+pdb_id.split("_")[0] + ".cumulative_acceptance.pdf")
                if not x_values:
                    x_values = make_figure.get_x_data(make_figure.labels[0])

                make_all_pdbids_figure.add_data(x_values, numpy.mean(make_figure.get_y_as_list(make_figure.labels), axis = 0), pdb_id)
                all_y_values.extend(make_figure.get_y_as_list(make_figure.labels))

            make_exp_figure.add_data(x_values, numpy.mean(all_y_values, axis=0), exp)

            make_all_pdbids_figure.fill_subplot("Cumulative MC Acceptance", make_all_pdbids_figure.labels, x_axis_label="Cycle", y_axis_label="Avg Cumulative Accepts")
            make_all_pdbids_figure.save_plot(self.plot_outdir+"/"+"all_pdbs_avg_cumulative_acceptance_"+exp+".pdf")

        make_exp_figure.fill_subplot("Cumulative MC Acceptance", make_exp_figure.labels, x_axis_label="Cycle", y_axis_label="Avg Cumulative Accepts", add_legend=True)
        make_exp_figure.save_plot(self.plot_outdir+"/"+"avg_cumulative_acceptance_"+"_".join([exp for exp in self.all_data])+".pdf")