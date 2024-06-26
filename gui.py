# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_RTA(object):
    def setupUi(self, RTA):
        RTA.setObjectName("RTA")
        RTA.resize(798, 568)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("images/logo_fih.jpg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        RTA.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(RTA)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(10, 0, 781, 371))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.implicit_gp = QtWidgets.QGroupBox(self.tab)
        self.implicit_gp.setGeometry(QtCore.QRect(20, 0, 731, 341))
        self.implicit_gp.setObjectName("implicit_gp")
        self.casematerial_gp = QtWidgets.QGroupBox(self.implicit_gp)
        self.casematerial_gp.setGeometry(QtCore.QRect(10, 19, 351, 101))
        self.casematerial_gp.setObjectName("casematerial_gp")
        self.btn_custom_case = QtWidgets.QPushButton(self.casematerial_gp)
        self.btn_custom_case.setGeometry(QtCore.QRect(190, 40, 121, 23))
        self.btn_custom_case.setObjectName("btn_custom_case")
        self.material_case = QtWidgets.QComboBox(self.casematerial_gp)
        self.material_case.setGeometry(QtCore.QRect(20, 40, 151, 22))
        self.material_case.setObjectName("material_case")
        self.material_case.addItem("")
        self.material_case.addItem("")
        self.material_case.addItem("")
        self.material_case.addItem("")
        self.material_case.addItem("")
        self.linermaterial_gp = QtWidgets.QGroupBox(self.implicit_gp)
        self.linermaterial_gp.setGeometry(QtCore.QRect(380, 20, 341, 101))
        self.linermaterial_gp.setObjectName("linermaterial_gp")
        self.btn_custom_insulator = QtWidgets.QPushButton(self.linermaterial_gp)
        self.btn_custom_insulator.setGeometry(QtCore.QRect(190, 40, 141, 23))
        self.btn_custom_insulator.setObjectName("btn_custom_insulator")
        self.material_insulator = QtWidgets.QComboBox(self.linermaterial_gp)
        self.material_insulator.setGeometry(QtCore.QRect(30, 40, 151, 22))
        self.material_insulator.setObjectName("material_insulator")
        self.material_insulator.addItem("")
        self.material_insulator.addItem("")
        self.material_insulator.addItem("")
        self.material_insulator.addItem("")
        self.material_insulator.addItem("")
        self.properties_implicit_gp = QtWidgets.QGroupBox(self.implicit_gp)
        self.properties_implicit_gp.setGeometry(QtCore.QRect(10, 130, 711, 201))
        self.properties_implicit_gp.setObjectName("properties_implicit_gp")
        self.label = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label.setGeometry(QtCore.QRect(20, 50, 121, 16))
        self.label.setObjectName("label")
        self.insulator_thk = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.insulator_thk.setGeometry(QtCore.QRect(20, 70, 141, 20))
        self.insulator_thk.setObjectName("insulator_thk")
        self.label_2 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_2.setGeometry(QtCore.QRect(20, 100, 121, 16))
        self.label_2.setObjectName("label_2")
        self.case_thk = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.case_thk.setGeometry(QtCore.QRect(20, 120, 141, 20))
        self.case_thk.setObjectName("case_thk")
        self.label_3 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_3.setGeometry(QtCore.QRect(20, 150, 131, 16))
        self.label_3.setObjectName("label_3")
        self.ri = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.ri.setGeometry(QtCore.QRect(20, 170, 141, 20))
        self.ri.setObjectName("ri")
        self.label_4 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_4.setGeometry(QtCore.QRect(170, 50, 121, 16))
        self.label_4.setObjectName("label_4")
        self.r_steps = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.r_steps.setGeometry(QtCore.QRect(170, 70, 151, 20))
        self.r_steps.setObjectName("r_steps")
        self.hm = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.hm.setGeometry(QtCore.QRect(170, 170, 151, 20))
        self.hm.setObjectName("hm")
        self.label_5 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_5.setGeometry(QtCore.QRect(170, 150, 141, 16))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_6.setGeometry(QtCore.QRect(330, 150, 141, 16))
        self.label_6.setObjectName("label_6")
        self.Ta = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.Ta.setGeometry(QtCore.QRect(330, 170, 151, 20))
        self.Ta.setObjectName("Ta")
        self.label_8 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_8.setGeometry(QtCore.QRect(330, 100, 141, 16))
        self.label_8.setObjectName("label_8")
        self.Tc = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.Tc.setGeometry(QtCore.QRect(330, 120, 151, 20))
        self.Tc.setObjectName("Tc")
        self.label_9 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_9.setGeometry(QtCore.QRect(170, 100, 141, 16))
        self.label_9.setObjectName("label_9")
        self.t_steps = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.t_steps.setGeometry(QtCore.QRect(170, 120, 151, 20))
        self.t_steps.setObjectName("t_steps")
        self.label_10 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_10.setGeometry(QtCore.QRect(330, 50, 141, 16))
        self.label_10.setObjectName("label_10")
        self.burn_time = QtWidgets.QLineEdit(self.properties_implicit_gp)
        self.burn_time.setGeometry(QtCore.QRect(331, 70, 151, 20))
        self.burn_time.setObjectName("burn_time")
        self.run_implicit = QtWidgets.QPushButton(self.properties_implicit_gp)
        self.run_implicit.setGeometry(QtCore.QRect(510, 100, 171, 51))
        self.run_implicit.setObjectName("run_implicit")
        self.label_11 = QtWidgets.QLabel(self.properties_implicit_gp)
        self.label_11.setGeometry(QtCore.QRect(120, 20, 41, 16))
        self.label_11.setObjectName("label_11")
        self.type_analysis_case = QtWidgets.QComboBox(self.properties_implicit_gp)
        self.type_analysis_case.setGeometry(QtCore.QRect(170, 20, 151, 22))
        self.type_analysis_case.setObjectName("type_analysis_case")
        self.type_analysis_case.addItem("")
        self.type_analysis_case.addItem("")
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.bulkhead_gp = QtWidgets.QGroupBox(self.tab_2)
        self.bulkhead_gp.setGeometry(QtCore.QRect(10, 10, 751, 311))
        self.bulkhead_gp.setObjectName("bulkhead_gp")
        self.bulkhead_material_gp = QtWidgets.QGroupBox(self.bulkhead_gp)
        self.bulkhead_material_gp.setGeometry(QtCore.QRect(10, 19, 351, 101))
        self.bulkhead_material_gp.setObjectName("bulkhead_material_gp")
        self.material_bulk = QtWidgets.QComboBox(self.bulkhead_material_gp)
        self.material_bulk.setGeometry(QtCore.QRect(20, 40, 151, 22))
        self.material_bulk.setObjectName("material_bulk")
        self.material_bulk.addItem("")
        self.material_bulk.addItem("")
        self.material_bulk.addItem("")
        self.material_bulk.addItem("")
        self.material_bulk.addItem("")
        self.btn_custom_case_2 = QtWidgets.QPushButton(self.bulkhead_material_gp)
        self.btn_custom_case_2.setGeometry(QtCore.QRect(210, 40, 121, 23))
        self.btn_custom_case_2.setObjectName("btn_custom_case_2")
        self.liner_material_gp = QtWidgets.QGroupBox(self.bulkhead_gp)
        self.liner_material_gp.setGeometry(QtCore.QRect(370, 20, 371, 101))
        self.liner_material_gp.setObjectName("liner_material_gp")
        self.btn_custom_insulator_2 = QtWidgets.QPushButton(self.liner_material_gp)
        self.btn_custom_insulator_2.setGeometry(QtCore.QRect(210, 40, 141, 23))
        self.btn_custom_insulator_2.setObjectName("btn_custom_insulator_2")
        self.material_insulator_2 = QtWidgets.QComboBox(self.liner_material_gp)
        self.material_insulator_2.setGeometry(QtCore.QRect(20, 40, 151, 22))
        self.material_insulator_2.setObjectName("material_insulator_2")
        self.material_insulator_2.addItem("")
        self.material_insulator_2.addItem("")
        self.material_insulator_2.addItem("")
        self.material_insulator_2.addItem("")
        self.material_insulator_2.addItem("")
        self.properties_bulk_gp = QtWidgets.QGroupBox(self.bulkhead_gp)
        self.properties_bulk_gp.setGeometry(QtCore.QRect(10, 120, 731, 181))
        self.properties_bulk_gp.setObjectName("properties_bulk_gp")
        self.label_29 = QtWidgets.QLabel(self.properties_bulk_gp)
        self.label_29.setGeometry(QtCore.QRect(10, 20, 121, 16))
        self.label_29.setObjectName("label_29")
        self.insulator_thk_4 = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.insulator_thk_4.setGeometry(QtCore.QRect(10, 40, 161, 20))
        self.insulator_thk_4.setObjectName("insulator_thk_4")
        self.label_30 = QtWidgets.QLabel(self.properties_bulk_gp)
        self.label_30.setGeometry(QtCore.QRect(10, 120, 121, 16))
        self.label_30.setObjectName("label_30")
        self.z = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.z.setGeometry(QtCore.QRect(10, 140, 161, 20))
        self.z.setObjectName("z")
        self.radius_inner = QtWidgets.QLabel(self.properties_bulk_gp)
        self.radius_inner.setGeometry(QtCore.QRect(180, 120, 131, 16))
        self.radius_inner.setObjectName("radius_inner")
        self.ri_4 = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.ri_4.setGeometry(QtCore.QRect(180, 140, 161, 20))
        self.ri_4.setTabletTracking(True)
        self.ri_4.setObjectName("ri_4")
        self.label_32 = QtWidgets.QLabel(self.properties_bulk_gp)
        self.label_32.setGeometry(QtCore.QRect(180, 20, 121, 16))
        self.label_32.setObjectName("label_32")
        self.r_steps_4 = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.r_steps_4.setGeometry(QtCore.QRect(180, 40, 161, 20))
        self.r_steps_4.setObjectName("r_steps_4")
        self.hm_3 = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.hm_3.setGeometry(QtCore.QRect(180, 90, 161, 20))
        self.hm_3.setObjectName("hm_3")
        self.label_33 = QtWidgets.QLabel(self.properties_bulk_gp)
        self.label_33.setGeometry(QtCore.QRect(180, 70, 141, 16))
        self.label_33.setObjectName("label_33")
        self.label_34 = QtWidgets.QLabel(self.properties_bulk_gp)
        self.label_34.setGeometry(QtCore.QRect(350, 120, 141, 16))
        self.label_34.setObjectName("label_34")
        self.Ta_3 = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.Ta_3.setGeometry(QtCore.QRect(350, 140, 161, 20))
        self.Ta_3.setObjectName("Ta_3")
        self.label_35 = QtWidgets.QLabel(self.properties_bulk_gp)
        self.label_35.setGeometry(QtCore.QRect(350, 70, 141, 16))
        self.label_35.setObjectName("label_35")
        self.Tc_3 = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.Tc_3.setGeometry(QtCore.QRect(350, 90, 161, 20))
        self.Tc_3.setObjectName("Tc_3")
        self.label_37 = QtWidgets.QLabel(self.properties_bulk_gp)
        self.label_37.setGeometry(QtCore.QRect(10, 70, 141, 16))
        self.label_37.setObjectName("label_37")
        self.burn_time_4 = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.burn_time_4.setGeometry(QtCore.QRect(10, 90, 161, 20))
        self.burn_time_4.setObjectName("burn_time_4")
        self.radius_outer = QtWidgets.QLineEdit(self.properties_bulk_gp)
        self.radius_outer.setGeometry(QtCore.QRect(350, 40, 161, 20))
        self.radius_outer.setObjectName("radius_outer")
        self.r_outer = QtWidgets.QLabel(self.properties_bulk_gp)
        self.r_outer.setGeometry(QtCore.QRect(350, 20, 131, 16))
        self.r_outer.setObjectName("r_outer")
        self.run_bulk = QtWidgets.QPushButton(self.properties_bulk_gp)
        self.run_bulk.setGeometry(QtCore.QRect(540, 70, 161, 61))
        self.run_bulk.setObjectName("run_bulk")
        self.tabWidget.addTab(self.tab_2, "")
        self.convective_tab = QtWidgets.QWidget()
        self.convective_tab.setObjectName("convective_tab")
        self.tau_logo_2 = QtWidgets.QLabel(self.convective_tab)
        self.tau_logo_2.setGeometry(QtCore.QRect(460, 0, 241, 111))
        self.tau_logo_2.setObjectName("tau_logo_2")
        self.label_25 = QtWidgets.QLabel(self.convective_tab)
        self.label_25.setGeometry(QtCore.QRect(20, 260, 341, 41))
        self.label_25.setObjectName("label_25")
        self.label_26 = QtWidgets.QLabel(self.convective_tab)
        self.label_26.setGeometry(QtCore.QRect(20, 290, 301, 20))
        self.label_26.setObjectName("label_26")
        self.label_27 = QtWidgets.QLabel(self.convective_tab)
        self.label_27.setGeometry(QtCore.QRect(420, 240, 341, 41))
        self.label_27.setObjectName("label_27")
        self.label_28 = QtWidgets.QLabel(self.convective_tab)
        self.label_28.setGeometry(QtCore.QRect(420, 270, 341, 16))
        self.label_28.setObjectName("label_28")
        self.label_31 = QtWidgets.QLabel(self.convective_tab)
        self.label_31.setGeometry(QtCore.QRect(420, 100, 341, 41))
        self.label_31.setObjectName("label_31")
        self.label_36 = QtWidgets.QLabel(self.convective_tab)
        self.label_36.setGeometry(QtCore.QRect(510, 120, 181, 16))
        self.label_36.setObjectName("label_36")
        self.label_38 = QtWidgets.QLabel(self.convective_tab)
        self.label_38.setGeometry(QtCore.QRect(510, 140, 181, 16))
        self.label_38.setObjectName("label_38")
        self.label_39 = QtWidgets.QLabel(self.convective_tab)
        self.label_39.setGeometry(QtCore.QRect(510, 160, 211, 16))
        self.label_39.setObjectName("label_39")
        self.label_40 = QtWidgets.QLabel(self.convective_tab)
        self.label_40.setGeometry(QtCore.QRect(510, 180, 211, 16))
        self.label_40.setObjectName("label_40")
        self.label_41 = QtWidgets.QLabel(self.convective_tab)
        self.label_41.setGeometry(QtCore.QRect(510, 200, 211, 16))
        self.label_41.setObjectName("label_41")
        self.label_42 = QtWidgets.QLabel(self.convective_tab)
        self.label_42.setGeometry(QtCore.QRect(510, 220, 211, 16))
        self.label_42.setObjectName("label_42")
        self.label_43 = QtWidgets.QLabel(self.convective_tab)
        self.label_43.setGeometry(QtCore.QRect(420, 290, 331, 16))
        self.label_43.setObjectName("label_43")
        self.label_44 = QtWidgets.QLabel(self.convective_tab)
        self.label_44.setGeometry(QtCore.QRect(420, 310, 331, 16))
        self.label_44.setObjectName("label_44")
        self.groupBox = QtWidgets.QGroupBox(self.convective_tab)
        self.groupBox.setGeometry(QtCore.QRect(10, 10, 361, 251))
        self.groupBox.setObjectName("groupBox")
        self.label_23 = QtWidgets.QLabel(self.groupBox)
        self.label_23.setGeometry(QtCore.QRect(190, 20, 111, 16))
        self.label_23.setObjectName("label_23")
        self.propellant_cp = QtWidgets.QLineEdit(self.groupBox)
        self.propellant_cp.setGeometry(QtCore.QRect(10, 50, 161, 20))
        self.propellant_cp.setObjectName("propellant_cp")
        self.label_21 = QtWidgets.QLabel(self.groupBox)
        self.label_21.setGeometry(QtCore.QRect(190, 90, 101, 16))
        self.label_21.setObjectName("label_21")
        self.hm_results = QtWidgets.QLineEdit(self.groupBox)
        self.hm_results.setGeometry(QtCore.QRect(190, 170, 161, 20))
        self.hm_results.setObjectName("hm_results")
        self.motor_length = QtWidgets.QLineEdit(self.groupBox)
        self.motor_length.setGeometry(QtCore.QRect(190, 120, 161, 20))
        self.motor_length.setObjectName("motor_length")
        self.label_18 = QtWidgets.QLabel(self.groupBox)
        self.label_18.setGeometry(QtCore.QRect(10, 80, 151, 31))
        self.label_18.setObjectName("label_18")
        self.label_22 = QtWidgets.QLabel(self.groupBox)
        self.label_22.setGeometry(QtCore.QRect(10, 150, 161, 16))
        self.label_22.setObjectName("label_22")
        self.hm_button = QtWidgets.QPushButton(self.groupBox)
        self.hm_button.setGeometry(QtCore.QRect(10, 210, 161, 23))
        self.hm_button.setObjectName("hm_button")
        self.burn_time_hm = QtWidgets.QLineEdit(self.groupBox)
        self.burn_time_hm.setGeometry(QtCore.QRect(10, 170, 161, 20))
        self.burn_time_hm.setObjectName("burn_time_hm")
        self.ri_hm = QtWidgets.QLineEdit(self.groupBox)
        self.ri_hm.setGeometry(QtCore.QRect(10, 120, 161, 20))
        self.ri_hm.setObjectName("ri_hm")
        self.hm_pdfreport = QtWidgets.QPushButton(self.groupBox)
        self.hm_pdfreport.setGeometry(QtCore.QRect(190, 210, 161, 23))
        self.hm_pdfreport.setObjectName("hm_pdfreport")
        self.label_24 = QtWidgets.QLabel(self.groupBox)
        self.label_24.setGeometry(QtCore.QRect(190, 150, 101, 16))
        self.label_24.setObjectName("label_24")
        self.label_7 = QtWidgets.QLabel(self.groupBox)
        self.label_7.setGeometry(QtCore.QRect(10, 20, 161, 16))
        self.label_7.setObjectName("label_7")
        self.propellant_mass = QtWidgets.QLineEdit(self.groupBox)
        self.propellant_mass.setGeometry(QtCore.QRect(190, 50, 161, 20))
        self.propellant_mass.setObjectName("propellant_mass")
        self.tabWidget.addTab(self.convective_tab, "")
        self.tau_logo = QtWidgets.QLabel(self.centralwidget)
        self.tau_logo.setGeometry(QtCore.QRect(250, 400, 311, 121))
        self.tau_logo.setObjectName("tau_logo")
        RTA.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(RTA)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 798, 21))
        self.menubar.setObjectName("menubar")
        self.menuAnalysis = QtWidgets.QMenu(self.menubar)
        self.menuAnalysis.setObjectName("menuAnalysis")
        self.menuResults = QtWidgets.QMenu(self.menubar)
        self.menuResults.setObjectName("menuResults")
        self.menuGenerate_output_file = QtWidgets.QMenu(self.menuResults)
        self.menuGenerate_output_file.setObjectName("menuGenerate_output_file")
        self.menuGenerate_graphs = QtWidgets.QMenu(self.menuResults)
        self.menuGenerate_graphs.setObjectName("menuGenerate_graphs")
        RTA.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(RTA)
        self.statusbar.setObjectName("statusbar")
        RTA.setStatusBar(self.statusbar)
        self.actionOpen = QtWidgets.QAction(RTA)
        self.actionOpen.setObjectName("actionOpen")
        self.actionSave = QtWidgets.QAction(RTA)
        self.actionSave.setObjectName("actionSave")
        self.generate_output_2 = QtWidgets.QAction(RTA)
        self.generate_output_2.setObjectName("generate_output_2")
        self.generate_output_3 = QtWidgets.QAction(RTA)
        self.generate_output_3.setObjectName("generate_output_3")
        self.generate_graphs_2 = QtWidgets.QAction(RTA)
        self.generate_graphs_2.setObjectName("generate_graphs_2")
        self.generate_graphs_3 = QtWidgets.QAction(RTA)
        self.generate_graphs_3.setObjectName("generate_graphs_3")
        self.menuAnalysis.addAction(self.actionOpen)
        self.menuAnalysis.addAction(self.actionSave)
        self.menuAnalysis.addSeparator()
        self.menuGenerate_output_file.addAction(self.generate_output_2)
        self.menuGenerate_output_file.addAction(self.generate_output_3)
        self.menuGenerate_graphs.addAction(self.generate_graphs_2)
        self.menuGenerate_graphs.addAction(self.generate_graphs_3)
        self.menuResults.addAction(self.menuGenerate_output_file.menuAction())
        self.menuResults.addSeparator()
        self.menuResults.addAction(self.menuGenerate_graphs.menuAction())
        self.menubar.addAction(self.menuAnalysis.menuAction())
        self.menubar.addAction(self.menuResults.menuAction())

        self.retranslateUi(RTA)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(RTA)
        RTA.setTabOrder(self.insulator_thk, self.case_thk)
        RTA.setTabOrder(self.case_thk, self.ri)
        RTA.setTabOrder(self.ri, self.r_steps)
        RTA.setTabOrder(self.r_steps, self.t_steps)
        RTA.setTabOrder(self.t_steps, self.hm)
        RTA.setTabOrder(self.hm, self.burn_time)
        RTA.setTabOrder(self.burn_time, self.Tc)
        RTA.setTabOrder(self.Tc, self.Ta)
        RTA.setTabOrder(self.Ta, self.run_implicit)
        RTA.setTabOrder(self.run_implicit, self.btn_custom_case)
        RTA.setTabOrder(self.btn_custom_case, self.btn_custom_insulator)
        RTA.setTabOrder(self.btn_custom_insulator, self.tabWidget)
        RTA.setTabOrder(self.tabWidget, self.insulator_thk_4)
        RTA.setTabOrder(self.insulator_thk_4, self.burn_time_4)
        RTA.setTabOrder(self.burn_time_4, self.z)
        RTA.setTabOrder(self.z, self.r_steps_4)
        RTA.setTabOrder(self.r_steps_4, self.hm_3)
        RTA.setTabOrder(self.hm_3, self.ri_4)
        RTA.setTabOrder(self.ri_4, self.radius_outer)
        RTA.setTabOrder(self.radius_outer, self.Tc_3)
        RTA.setTabOrder(self.Tc_3, self.Ta_3)
        RTA.setTabOrder(self.Ta_3, self.run_bulk)
        RTA.setTabOrder(self.run_bulk, self.type_analysis_case)
        RTA.setTabOrder(self.type_analysis_case, self.material_case)
        RTA.setTabOrder(self.material_case, self.material_insulator)
        RTA.setTabOrder(self.material_insulator, self.material_bulk)
        RTA.setTabOrder(self.material_bulk, self.btn_custom_case_2)
        RTA.setTabOrder(self.btn_custom_case_2, self.btn_custom_insulator_2)
        RTA.setTabOrder(self.btn_custom_insulator_2, self.material_insulator_2)
        RTA.setTabOrder(self.material_insulator_2, self.propellant_cp)
        RTA.setTabOrder(self.propellant_cp, self.hm_results)
        RTA.setTabOrder(self.hm_results, self.motor_length)
        RTA.setTabOrder(self.motor_length, self.hm_button)
        RTA.setTabOrder(self.hm_button, self.burn_time_hm)
        RTA.setTabOrder(self.burn_time_hm, self.ri_hm)
        RTA.setTabOrder(self.ri_hm, self.hm_pdfreport)
        RTA.setTabOrder(self.hm_pdfreport, self.propellant_mass)

    def retranslateUi(self, RTA):
        _translate = QtCore.QCoreApplication.translate
        RTA.setWindowTitle(_translate("RTA", "Rocket Thermal Analysis"))
        self.implicit_gp.setTitle(_translate("RTA", "Inputs"))
        self.casematerial_gp.setTitle(_translate("RTA", "Case Materials"))
        self.btn_custom_case.setText(_translate("RTA", "Custom case material"))
        self.material_case.setItemText(0, _translate("RTA", "Aluminium 6061-T6"))
        self.material_case.setItemText(1, _translate("RTA", "Stainless Steel 304"))
        self.material_case.setItemText(2, _translate("RTA", "Stainless Steel 316"))
        self.material_case.setItemText(3, _translate("RTA", "Steel 1010"))
        self.material_case.setItemText(4, _translate("RTA", "Carbon Fiber"))
        self.linermaterial_gp.setTitle(_translate("RTA", "Insulator materials"))
        self.btn_custom_insulator.setText(_translate("RTA", "Custom insulation material"))
        self.material_insulator.setItemText(0, _translate("RTA", "EPDM"))
        self.material_insulator.setItemText(1, _translate("RTA", "NBR"))
        self.material_insulator.setItemText(2, _translate("RTA", "Paper"))
        self.material_insulator.setItemText(3, _translate("RTA", "FiberGlass"))
        self.material_insulator.setItemText(4, _translate("RTA", "Phenolic Paper"))
        self.properties_implicit_gp.setTitle(_translate("RTA", "Properties"))
        self.label.setText(_translate("RTA", "Insulator Thickness [m]:"))
        self.label_2.setText(_translate("RTA", "Case Thickness [m]:"))
        self.label_3.setText(_translate("RTA", "Inner cylinder radius [m]:"))
        self.label_4.setText(_translate("RTA", "Radial divisions:"))
        self.label_5.setText(_translate("RTA", "Convection coeff [W/m2 -K]:"))
        self.label_6.setText(_translate("RTA", "Initial Temperature [K]:"))
        self.label_8.setText(_translate("RTA", "Combustion Temperature [K]:"))
        self.label_9.setText(_translate("RTA", "Time steps:"))
        self.label_10.setText(_translate("RTA", "Burn time [s]:"))
        self.run_implicit.setText(_translate("RTA", "Run analysis"))
        self.label_11.setText(_translate("RTA", "Solver:"))
        self.type_analysis_case.setItemText(0, _translate("RTA", "Implicit"))
        self.type_analysis_case.setItemText(1, _translate("RTA", "Explicit"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("RTA", "Motor case"))
        self.bulkhead_gp.setTitle(_translate("RTA", "Inputs [Explicit Method]:"))
        self.bulkhead_material_gp.setTitle(_translate("RTA", "Bulkhead Materials"))
        self.material_bulk.setItemText(0, _translate("RTA", "Aluminium 6061-T6"))
        self.material_bulk.setItemText(1, _translate("RTA", "Stainless Steel 304"))
        self.material_bulk.setItemText(2, _translate("RTA", "Stainless Steel 316"))
        self.material_bulk.setItemText(3, _translate("RTA", "Steel 1010"))
        self.material_bulk.setItemText(4, _translate("RTA", "Carbon Fiber"))
        self.btn_custom_case_2.setText(_translate("RTA", "Custom case material"))
        self.liner_material_gp.setTitle(_translate("RTA", "Insulator materials"))
        self.btn_custom_insulator_2.setText(_translate("RTA", "Custom insulation material"))
        self.material_insulator_2.setItemText(0, _translate("RTA", "EPDM"))
        self.material_insulator_2.setItemText(1, _translate("RTA", "NBR"))
        self.material_insulator_2.setItemText(2, _translate("RTA", "Paper"))
        self.material_insulator_2.setItemText(3, _translate("RTA", "FiberGlass"))
        self.material_insulator_2.setItemText(4, _translate("RTA", "Phenolic Paper"))
        self.properties_bulk_gp.setTitle(_translate("RTA", "Properties"))
        self.label_29.setText(_translate("RTA", "Insulator Thickness [m]:"))
        self.label_30.setText(_translate("RTA", "Bulkhead Thickness [m]:"))
        self.radius_inner.setText(_translate("RTA", "Chamber inner radius [m]:"))
        self.label_32.setText(_translate("RTA", "Radial Sections:"))
        self.label_33.setText(_translate("RTA", "Convection coeff [W/m2 -K]:"))
        self.label_34.setText(_translate("RTA", "Initial Temperature [K]:"))
        self.label_35.setText(_translate("RTA", "Combustion Temperature [K]:"))
        self.label_37.setText(_translate("RTA", "Burn time:"))
        self.r_outer.setText(_translate("RTA", "Chamber outer radius [m]:"))
        self.run_bulk.setText(_translate("RTA", "Run analysis"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("RTA", "Bulkhead"))
        self.tau_logo_2.setText(_translate("RTA", "<html><head/><body><p><img src=\":/images/equations.PNG\" width=\"250\" height=\"100\"/></p></body></html>"))
        self.label_25.setText(_translate("RTA", "The convective heat transfer coefficient was calculated using the "))
        self.label_26.setText(_translate("RTA", "formula shown in right hand side equation:"))
        self.label_27.setText(_translate("RTA", "The specific heat of common propellants:"))
        self.label_28.setText(_translate("RTA", "KNSU 65/35: 1690 J/Kg-K                      KNSB 65/35: 1740 J/Kg-K"))
        self.label_31.setText(_translate("RTA", "Where,"))
        self.label_36.setText(_translate("RTA", "Cp -  Propellant specific heat"))
        self.label_38.setText(_translate("RTA", "G   -  Average mass flux"))
        self.label_39.setText(_translate("RTA", "Di  -  Inner diameter of combution chamber"))
        self.label_40.setText(_translate("RTA", "L   -  Length of the motor"))
        self.label_41.setText(_translate("RTA", "t   -  Burn time"))
        self.label_42.setText(_translate("RTA", "mp - Propellant mass"))
        self.label_43.setText(_translate("RTA", "KNDX 65/35: 1700 J/Kg-K"))
        self.label_44.setText(_translate("RTA", "See more in: Mark’s Standard Handbook for Mechanical Engineers."))
        self.groupBox.setTitle(_translate("RTA", "Inputs needed"))
        self.label_23.setText(_translate("RTA", "Propellant mass [Kg]:"))
        self.label_21.setText(_translate("RTA", "Motor Length [m]:"))
        self.label_18.setText(_translate("RTA", "Inner combustion \n"
"chamber diameter [m]:"))
        self.label_22.setText(_translate("RTA", "Burn Time [s]:"))
        self.hm_button.setText(_translate("RTA", "Calculate"))
        self.hm_pdfreport.setText(_translate("RTA", "Generate PDF report"))
        self.label_24.setText(_translate("RTA", "Result [W/m^2-K]:"))
        self.label_7.setText(_translate("RTA", "Propellant Specific Heat [J/Kg-K]:"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.convective_tab), _translate("RTA", "Convective Coefficient"))
        self.tau_logo.setText(_translate("RTA", "<html><head/><body><p><img src=\":/images/logo_horizontal.png\"  width=\"300\" height=\"95\" /></p></body></html>"))
        self.menuAnalysis.setTitle(_translate("RTA", "Analysis"))
        self.menuResults.setTitle(_translate("RTA", "Results"))
        self.menuGenerate_output_file.setTitle(_translate("RTA", "Generate output file"))
        self.menuGenerate_graphs.setTitle(_translate("RTA", "Generate graphs"))
        self.actionOpen.setText(_translate("RTA", "Open"))
        self.actionSave.setText(_translate("RTA", "Save"))
        self.generate_output_2.setText(_translate("RTA", "From Motor Case Analysis"))
        self.generate_output_3.setText(_translate("RTA", "From Bulkhead Analysis"))
        self.generate_graphs_2.setText(_translate("RTA", "From Motor Case Analysis"))
        self.generate_graphs_3.setText(_translate("RTA", "From Bulkhead Analysis"))
from images import graphical_rc
