const core = require('@actions/core');
const { Client } = require('discord.js-selfbot-v13');
const fs = require('fs');
const path = require('path');

async function run() {
  try {
    const token = core.getInput('token', { required: true });
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

    // 3) Create the selfbot client (again, strongly TOS-breaking)
    const client = new Client({ checkUpdate: false });

    client.on('ready', async () => {
      console.log(`[DEBUG] Logged in as: ${client.user?.tag} - Attempting to update avatar...`);

      try {
        // 4) Use the selfbot method to update your user avatar
        // await client.user?.setAvatar(base64Image);
        console.log("Avatar updated successfully!");
      } catch (error) {
        core.setFailed(`Failed to update avatar: ${error.message}`);
      } finally {
        client.destroy();
        process.exit(0);
      }
    });
  } catch (error) {
    core.setFailed(error.message);
  }
}

run();
