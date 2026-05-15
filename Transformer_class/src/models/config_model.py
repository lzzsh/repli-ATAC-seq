from transformers import PretrainedConfig


class RepliformerConfig(PretrainedConfig):
    model_type = "repliformer"

    def __init__(
        self,
        # ── Conv trunk ────────────────────────────────────────────────────────
        stem_channels: int = 768,
        stem_kernel_size: int = 15,
        tower_filter_list: list = None,   # None → default Enformer schedule
        bn_momentum: float = 0.1,
        # ── Pooling ───────────────────────────────────────────────────────────
        pool_type: str = "attn",        # "attn" | "avg_late" (avg for last 3 tower pools)
        # ── Transformer tower ─────────────────────────────────────────────────
        d_model: int = 1536,
        n_heads: int = 8,
        dim_key: int = 64,
        n_layers: int = 11,
        dropout: float = 0.4,
        attn_dropout: float = 0.05,
        pos_dropout: float = 0.01,
        # ── Final pointwise head ──────────────────────────────────────────────
        pointwise_channels: int = 3072,
        final_dropout: float = 0.05,
        # ── Sequence / output geometry ────────────────────────────────────────
        crop_size: int = 320,           # bins cropped from each side before head
        n_output_bins: int = 3,         # RT classes per bin (ES, MS, LS)
        **kwargs,
    ):
        # default Enformer conv-tower filter schedule
        if tower_filter_list is None:
            tower_filter_list = [768, 768, 896, 1024, 1152, 1280, 1536]

        self.stem_channels = stem_channels
        self.stem_kernel_size = stem_kernel_size
        self.tower_filter_list = tower_filter_list
        self.pool_type = pool_type
        self.bn_momentum = bn_momentum

        self.d_model = d_model
        self.n_heads = n_heads
        self.dim_key = dim_key
        self.n_layers = n_layers
        self.dropout = dropout
        self.attn_dropout = attn_dropout
        self.pos_dropout = pos_dropout

        self.pointwise_channels = pointwise_channels
        self.final_dropout = final_dropout

        self.crop_size = crop_size
        self.n_output_bins = n_output_bins

        super().__init__(**kwargs)
