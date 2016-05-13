package workshop;

import java.awt.Button;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.Label;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import jv.number.PuDouble;
import jv.object.PsConfig;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jvx.project.PjWorkshop_IP;

public class MyWorkshop_IP extends PjWorkshop_IP implements ActionListener {

	protected Button m_bMakeRandomElementColors;
	protected Button m_bMakeRandomVertexColors;
	protected PuDouble m_xOff;

	protected Button m_calculateGenusButton;
	protected Label m_genusLabel;

	protected Button m_volumeButton;
	protected Label m_volumeLabel;

	protected Button m_shapeRegButton;
	protected Label m_minShapeRegLabel;
	protected Label m_maxShapeRegLabel;
	protected Label m_meanShapeRegLabel;
	protected Label m_devShapeRegLabel;

	protected Button m_valenceButton;
	protected Label m_minValenceLabel;
	protected Label m_maxValenceLabel;
	protected Label m_meanValenceLabel;
	protected Label m_devValenceLabel;
	
	MyWorkshop m_ws;
	
	public MyWorkshop_IP() {
		super();
		if(getClass() == MyWorkshop_IP.class)
			init();
	}
	
	public void init() {
		super.init();
		setTitle("My Workshop");
	}
	
	public String getNotice() {
		return "This text should explain what the workshop is about and how to use it.";
	}
	
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_ws = (MyWorkshop)parent;
	
		addSubTitle("Example of a subtitle");
		
		m_bMakeRandomElementColors = new Button("Random Element Colors");
		m_bMakeRandomElementColors.addActionListener(this);
		m_bMakeRandomVertexColors = new Button("Random Vertex Colors");
		m_bMakeRandomVertexColors.addActionListener(this);
		Panel panel1 = new Panel(new FlowLayout(FlowLayout.CENTER));
		panel1.add(m_bMakeRandomElementColors);
		panel1.add(m_bMakeRandomVertexColors);
		add(panel1);
		
		m_xOff = new PuDouble("X Offset");
		m_xOff.setDefBounds(-10,10,0.1,1);
		m_xOff.addUpdateListener(this);
		m_xOff.init();
		add(m_xOff.getInfoPanel());

		m_calculateGenusButton = new Button("Calculate Genus");
		m_calculateGenusButton.addActionListener(this);
		add(m_calculateGenusButton, FlowLayout.CENTER);
		m_genusLabel = new Label("Genus: ?");
		add(m_genusLabel, FlowLayout.CENTER);

		m_volumeButton = new Button("Calculate Volume");
		m_volumeButton.addActionListener(this);
		add(m_volumeButton, FlowLayout.CENTER);
		m_volumeLabel = new Label("Volume: ?");
		add(m_volumeLabel, FlowLayout.CENTER);

		m_shapeRegButton = new Button("Calculate Shape Regularity");
		m_shapeRegButton.addActionListener(this);
		add(m_shapeRegButton, FlowLayout.CENTER);
		m_minShapeRegLabel = new Label("Min shape regularity: ?");
		add(m_minShapeRegLabel, FlowLayout.CENTER);
		m_maxShapeRegLabel = new Label("Max shape regularity: ?");
		add(m_maxShapeRegLabel, FlowLayout.CENTER);
		m_meanShapeRegLabel = new Label("Average shape regularity: ?");
		add(m_meanShapeRegLabel, FlowLayout.CENTER);
		m_devShapeRegLabel = new Label("Standard deviation shape regularity: ?");
		add(m_devShapeRegLabel, FlowLayout.CENTER);

		m_valenceButton = new Button("Calculate Valence");
		m_valenceButton.addActionListener(this);
		add(m_valenceButton, FlowLayout.CENTER);
		m_minValenceLabel = new Label("Min valence: ?");
		add(m_minValenceLabel, FlowLayout.CENTER);
		m_maxValenceLabel = new Label("Max valence: ?");
		add(m_maxValenceLabel, FlowLayout.CENTER);
		m_meanValenceLabel = new Label("Average valence: ?");
		add(m_meanValenceLabel, FlowLayout.CENTER);
		m_devValenceLabel = new Label("Standard deviation valence: ?");
		add(m_devValenceLabel, FlowLayout.CENTER);
		
		validate();
	}
	
	
	public boolean update(Object event) {
		if (event == m_xOff) {
			m_ws.setXOff(m_xOff.getValue());
			m_ws.m_geom.update(m_ws.m_geom);
			return true;
		} else
			return super.update(event);
	}
	
	/**
	 * Handle action events fired by buttons etc.
	 */
	public void actionPerformed(ActionEvent event) {
		Object source = event.getSource();
		if (source == m_bMakeRandomElementColors) {
			m_ws.makeRandomElementColors();
			m_ws.m_geom.update(m_ws.m_geom);
			return;
		}
		else if (source == m_bMakeRandomVertexColors) {
			m_ws.makeRandomVertexColors();
			m_ws.m_geom.update(m_ws.m_geom);
			return;
		}
		else if (source == m_calculateGenusButton)
		{
			int genus = m_ws.calcGenus();
			m_genusLabel.setText("Genus: " + genus);
		}
		else if (source == m_volumeButton)
		{
			double volume = m_ws.calcVolume();
			m_volumeLabel.setText("Volume: " + volume);
		}
		else if (source == m_shapeRegButton)
		{
			double[] stats = m_ws.calcShapeReg();
			m_ws.makeRegularityElementColors();
			m_ws.m_geom.update(m_ws.m_geom);

			m_minShapeRegLabel.setText("Min shape regularity: " + stats[0]);
			m_maxShapeRegLabel.setText("Max shape regularity: " + stats[1]);
			m_meanShapeRegLabel.setText("Average shape regularity: " + stats[2]);
			m_devShapeRegLabel.setText("Standard deviation shape regularity: " + stats[3]);
		}
		else if (source == m_valenceButton)
		{
			double[] stats = m_ws.calcValence();
			m_minValenceLabel.setText("Min valence: " + stats[0]);
			m_maxValenceLabel.setText("Max valence: " + stats[1]);
			m_meanValenceLabel.setText("Average valence: " + stats[2]);
			m_devValenceLabel.setText("Standard deviation valence: " + stats[3]);
		}
	}
	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}
