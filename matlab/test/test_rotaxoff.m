function test_rotaxoff( fig, h, txt_pos)

    uicontrol(fig, 'Style','text',...
    'Position',txt_pos,...
    'String', sprintf('POS:%u', h.DataSource.Controls.CurrentFrame));

end

