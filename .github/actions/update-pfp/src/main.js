const core = require('@actions/core');
const { Client } = require('discord.js-selfbot-v13');
const fs = require('fs');
const path = require('path');
const captcha = require('2captcha');

async function run() {
  try {
    const token = process.env.DISCORD_TOKEN;
    const imageFile = core.getInput('imageFile', { required: true });

    console.log(`[DEBUG] Using token length: ${token.length} (masked)`);
    console.log(`[DEBUG] Requested imageFile: ${imageFile}`);

    const fullPath = path.join(process.env.GITHUB_WORKSPACE || '', imageFile);
    if (!fs.existsSync(fullPath)) {
      throw new Error(`Image file not found at: ${fullPath}`);
    }

    const fileData = fs.readFileSync(fullPath);
    const base64Image = `data:image/png;base64,${fileData.toString('base64')}`;

    console.log(`[DEBUG] Successfully read image file. Size: ${fileData.length} bytes`);

    const solver = new captcha.Solver(process.env.API_TOKEN_2CAPTCHA);

    const client = new Client({
      captchaSolver: function(captcha, UA) {
        return solver
            .hcaptcha(
                captcha.captcha_sitekey,
                'discord.com',
                {
                  invisible: 1,
                  userAgent: UA,
                  data: captcha.captcha_rqdata
                }
            ).then(res => res.data);
      },
    });

    client.on('ready', async () => {
      console.log(`[DEBUG] Logged in as: ${client.user?.tag} - Attempting to update avatar...`);

      try {
        await client.user?.setAvatar(base64Image);
        console.log("Avatar updated successfully!");
      } catch (error) {
        core.setFailed(`Failed to update avatar: ${error.message}`);
      } finally {
        client.destroy();
      }
    });

    await client.login(token);
  } catch (error) {
    core.setFailed(error.message);
  }
}

run();
